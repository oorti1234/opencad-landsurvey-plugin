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
    fn rejects_invalid_tolerance() {
        assert!(ambiguous_junctions(&[],0.).is_err());
    }
}
