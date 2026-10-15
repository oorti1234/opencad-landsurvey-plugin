//! Assemble independent isoline segments into ordered chains.
//! This module is experimental; CAD/host integration is deliberately separate.
use crate::contours::ContourSegment;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct ContourChain {
    pub elevation: f64,
    pub vertices: Vec<[f64; 2]>,
    pub closed: bool,
}

fn key(p: [f64; 2], tol: f64) -> Option<(i64, i64)> {
    if p.iter().any(|v| !v.is_finite()) { return None; }
    let x = (p[0] / tol).round();
    let y = (p[1] / tol).round();
    if x < i64::MIN as f64 || x > i64::MAX as f64 ||
       y < i64::MIN as f64 || y > i64::MAX as f64 { return None; }
    Some((x as i64, y as i64))
}

/// Join raw segments by elevation and endpoint coordinates. Returns an error
/// rather than guessing connectivity at branching or ambiguous nodes.
/// The caller chooses an XY tolerance expressed in drawing units.
pub fn stitch_segments(
    segments: &[ContourSegment],
    tolerance: f64,
) -> Result<Vec<ContourChain>, String> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err("tolerance must be finite and positive".into());
    }
    type K = (i64,i64);
    let mut by_level: BTreeMap<u64, Vec<(K,K)>> = BTreeMap::new();
    let mut coordinates: BTreeMap<(u64,K), [f64;2]> = BTreeMap::new();
    for seg in segments {
        if !seg.elevation.is_finite() { return Err("nonfinite elevation".into()); }
        let a = key(seg.start,tolerance).ok_or("invalid start point")?;
        let b = key(seg.end,tolerance).ok_or("invalid end point")?;
        if a == b { continue; }
        let lev = seg.elevation.to_bits();
        coordinates.entry((lev,a)).or_insert(seg.start);
        coordinates.entry((lev,b)).or_insert(seg.end);
        by_level.entry(lev).or_default().push((a,b));
    }
    let mut output = Vec::new();
    for (lev, edges) in by_level {
        let mut adj: BTreeMap<K, Vec<(usize,K)>> = BTreeMap::new();
        let mut unique = BTreeSet::new();
        let mut kept = Vec::new();
        for (a,b) in edges {
            let canon = if a < b {(a,b)} else {(b,a)};
            if !unique.insert(canon) { continue; }
            let id = kept.len();
            kept.push((a,b));
            adj.entry(a).or_default().push((id,b));
            adj.entry(b).or_default().push((id,a));
        }
        if let Some((node, neighbors)) = adj.iter().find(|(_, neighbors)| neighbors.len() > 2) {
            let xy = coordinates[&(lev, *node)];
            return Err(format!(
                "ambiguous contour junction (degree {}) at elevation {:.6}, X {:.6}, Y {:.6}; check exact-level TIN vertices or snapping tolerance",
                neighbors.len(), f64::from_bits(lev), xy[0], xy[1]
            ));
        }
        let mut visited = vec![false; kept.len()];
        // Open chains first, then cycles. Sorted map makes output deterministic.
        let starts: Vec<K> = adj.iter().filter(|(_,n)|n.len()==1).map(|(k,_)|*k)
            .chain(adj.keys().copied()).collect();
        for start in starts {
            if !adj[&start].iter().any(|(id,_)| !visited[*id]) { continue; }
            let mut path = vec![start];
            let mut curr = start;
            loop {
                let next = adj[&curr].iter().filter(|(id,_)|!visited[*id])
                    .min_by_key(|(_,other)|*other).copied();
                let Some((id,other)) = next else {break};
                visited[id]=true;
                path.push(other);
                curr=other;
                if curr==start {break}
            }
            let closed = path.len()>2 && path.first()==path.last();
            output.push(ContourChain {
                elevation: f64::from_bits(lev),
                vertices: path.iter().map(|k|coordinates[&(lev,*k)]).collect(),
                closed,
            });
        }
    }
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn segment(a:[f64;2],b:[f64;2])->ContourSegment{
        ContourSegment{elevation:5.,start:a,end:b}
    }
    #[test]
    fn joins_reversed_segments_into_open_chain(){
        let x=stitch_segments(&[
            segment([1.,0.],[0.,0.]),
            segment([2.,0.],[1.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert_eq!(x[0].vertices.len(),3);
        assert!(!x[0].closed);
    }
    #[test]
    fn identifies_closed_loop(){
        let x=stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([1.,0.],[1.,1.]),
            segment([1.,1.],[0.,1.]),
            segment([0.,1.],[0.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert!(x[0].closed);
        assert_eq!(x[0].vertices.len(),5);
    }
    #[test]
    fn ignores_reversed_duplicate_edge(){
        let x=stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([1.,0.],[0.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert_eq!(x[0].vertices.len(),2);
    }
    #[test]
    fn refuses_branch_points(){
        assert!(stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([0.,0.],[0.,1.]),
            segment([0.,0.],[-1.,0.])
        ],1e-8).is_err());
    }
    #[test]
    fn branch_diagnostic_identifies_level_and_position(){
        let err = stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([0.,0.],[0.,1.]),
            segment([0.,0.],[-1.,0.])
        ],1e-8).unwrap_err();
        assert!(err.contains("degree 3"));
        assert!(err.contains("elevation 5.000000"));
        assert!(err.contains("X 0.000000"));
    }
    #[test]
    fn rejects_invalid_tolerance(){
        assert!(stitch_segments(&[],0.).is_err());
    }
}
