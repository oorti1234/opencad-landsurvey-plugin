//! Non-destructive diagnostic of raw contour segment connectivity.
//! Identifies junctions where a contour cannot be paired unambiguously.
use crate::contours::ContourSegment;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct Junction {
    pub elevation: f64,
    pub xy: [f64; 2],
    pub degree: usize,
}

/// Analyze one or more levels without altering contour segments.
/// Duplicate/reversed edges count once, matching the stitching algorithm.
pub fn ambiguous_junctions(
    segments: &[ContourSegment],
    tolerance: f64,
) -> Result<Vec<Junction>, String> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err("invalid XY tolerance".into());
    }
    type Key = (u64, i64, i64);
    let mut coordinates: BTreeMap<Key, [f64; 2]> = BTreeMap::new();
    let mut edges: BTreeSet<(Key, Key)> = BTreeSet::new();
    for s in segments {
        if !s.elevation.is_finite() { return Err("invalid elevation".into()); }
        let make_key = |xy: [f64; 2]| -> Result<Key, String> {
            if xy.iter().any(|n| !n.is_finite()) { return Err("invalid XY".into()); }
            let x = (xy[0]/tolerance).round();
            let y = (xy[1]/tolerance).round();
            if x < i64::MIN as f64 || x >= i64::MAX as f64
                || y < i64::MIN as f64 || y >= i64::MAX as f64 {
                return Err("coordinate outside supported tolerance grid".into());
            }
            Ok((s.elevation.to_bits(), x as i64, y as i64))
        };
        let a = make_key(s.start)?;
        let b = make_key(s.end)?;
        if a == b { continue; }
        coordinates.entry(a).or_insert(s.start);
        coordinates.entry(b).or_insert(s.end);
        edges.insert(if a < b {(a,b)} else {(b,a)});
    }
    let mut degrees: BTreeMap<Key, usize> = BTreeMap::new();
    for (a,b) in edges {
        *degrees.entry(a).or_default() += 1;
        *degrees.entry(b).or_default() += 1;
    }
    Ok(degrees.into_iter().filter_map(|(k,degree)| {
        (degree > 2).then(|| Junction {
            elevation: f64::from_bits(k.0),
            xy: coordinates[&k],
            degree,
        })
    }).collect())
}

/// Classify a level-crossing node as a discrete saddle only if its
/// surrounding triangle fan is complete and adjacent elevations alternate.
/// Returns false for boundaries, flat neighboring vertices, or incomplete fans.
pub fn is_exact_level_saddle(
    surface: &crate::surface::Surface,
    xy: [f64; 2],
    elevation: f64,
    tolerance: f64,
) -> bool {
    if !tolerance.is_finite() || tolerance <= 0.0 || !elevation.is_finite()
        || xy.iter().any(|v| !v.is_finite()) { return false; }
    let candidates: Vec<usize> = surface.nodes.iter().enumerate()
        .filter(|(_,p)| (p[0]-xy[0]).hypot(p[1]-xy[1]) <= tolerance
            && (p[2]-elevation).abs() <= 1e-9)
        .map(|(i,_)|i).collect();
    if candidates.len()!=1 { return false; }
    let center=candidates[0];
    let mut neighbors = BTreeSet::new();
    let mut fan_edges = BTreeSet::new();
    let mut incident = 0usize;
    for tri in &surface.triangles {
        if tri.iter().any(|&i|i>=surface.nodes.len()) {return false;}
        if !tri.contains(&center) {continue;}
        let other: Vec<usize> = tri.iter().copied().filter(|&i|i!=center).collect();
        if other.len()!=2 || other[0]==other[1] {return false;}
        neighbors.extend(other.iter().copied());
        fan_edges.insert(if other[0]<other[1] {(other[0],other[1])}
            else {(other[1],other[0])});
        incident += 1;
    }
    if neighbors.len()<4 || incident!=neighbors.len() || fan_edges.len()!=neighbors.len() {
        return false;
    }
    let p=surface.nodes[center];
    let mut ordered: Vec<(f64,i8,usize)> = Vec::new();
    for n in neighbors {
        let q=surface.nodes[n];
        if q.iter().any(|x|!x.is_finite()) || (q[2]-elevation).abs()<=1e-9 {return false;}
        let dx=q[0]-p[0]; let dy=q[1]-p[1];
        if dx.hypot(dy)<=tolerance {return false;}
        ordered.push((dy.atan2(dx),if q[2]>elevation {1} else {-1},n));
    }
    ordered.sort_by(|a,b|a.0.total_cmp(&b.0));
    let mut flips=0usize;
    for i in 0..ordered.len() {
        let a=ordered[i].2;
        let b=ordered[(i+1)%ordered.len()].2;
        let edge=if a<b {(a,b)} else {(b,a)};
        if !fan_edges.contains(&edge) {return false;}
        if ordered[i].1!=ordered[(i+1)%ordered.len()].1 {flips+=1;}
    }
    flips>=4
}

#[cfg(test)]
mod tests {
    use super::*;
    fn seg(a:[f64;2],b:[f64;2])->ContourSegment {
        ContourSegment{elevation:41.5,start:a,end:b}
    }
    #[test]
    fn pinpoints_degree_four_at_exact_level() {
        let center = [355420.383,928758.559];
        let segments = [
            seg(center,[355421.,928758.559]),
            seg(center,[355419.,928758.559]),
            seg(center,[355420.383,928759.]),
            seg(center,[355420.383,928757.]),
        ];
        let problems = ambiguous_junctions(&segments,1e-6).unwrap();
        assert_eq!(problems.len(),1);
        assert_eq!(problems[0].degree,4);
        assert_eq!(problems[0].elevation,41.5);
        assert_eq!(problems[0].xy,center);
    }
    #[test]
    fn reversed_duplicate_is_not_new_branch() {
        let problems = ambiguous_junctions(&[
            seg([0.,0.],[1.,0.]), seg([1.,0.],[0.,0.]),
            seg([0.,0.],[0.,1.]),
        ],1e-6).unwrap();
        assert!(problems.is_empty());
    }
    #[test]
    fn detects_alternating_saddle_in_closed_fan() {
        let surface=crate::surface::Surface {
            nodes:vec![[0.,0.,41.5],[1.,0.,42.],[0.,1.,41.],[-1.,0.,42.],[0.,-1.,41.]],
            triangles:vec![[0,1,2],[0,2,3],[0,3,4],[0,4,1]],
        };
        assert!(is_exact_level_saddle(&surface,[0.,0.],41.5,1e-6));
    }
    #[test]
    fn boundary_vertex_is_not_classified_as_saddle() {
        let surface=crate::surface::Surface {
            nodes:vec![[0.,0.,41.5],[1.,0.,42.],[0.,1.,41.],[-1.,0.,42.],[0.,-1.,41.]],
            triangles:vec![[0,1,2],[0,2,3],[0,3,4]],
        };
        assert!(!is_exact_level_saddle(&surface,[0.,0.],41.5,1e-6));
    }
    #[test]
    fn rejects_invalid_tolerance() {
        assert!(ambiguous_junctions(&[],0.).is_err());
    }
}
