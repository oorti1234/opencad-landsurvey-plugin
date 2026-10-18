//! Experimental, host-independent TIN isoline segment generation.
//!
//! This is an initial building block, not yet a stitched-polyline or CAD command.
//! Coordinates are [Easting, Northing, Elevation], in the drawing's units.
//! Coplanar (flat-at-level) triangles are intentionally omitted for now;
//! their polygon boundaries require separate topology handling.

use crate::surface::Surface;

/// One contour segment: both endpoints lie on the same elevation level.
#[derive(Debug, Clone, PartialEq)]
pub struct ContourSegment {
    pub elevation: f64,
    pub start: [f64; 2],
    pub end: [f64; 2],
}

/// Generate raw isoline segments from a TIN at one level.
///
/// Uses a consistent half-open intersection convention: edges include their
/// low endpoint but exclude their high endpoint. Exact-level vertex ties are
/// handled by parity classification (< level versus >= level), suppressing
/// duplicate vertex crossings within a triangle. Flat edges are not emitted.
pub fn segments_at_level(surface: &Surface, level: f64) -> Vec<ContourSegment> {
    if !level.is_finite() { return Vec::new(); }
    let mut result = Vec::new();
    for &tri in &surface.triangles {
        if tri.iter().any(|&i| i >= surface.nodes.len()) { continue; }
        let p = [surface.nodes[tri[0]], surface.nodes[tri[1]], surface.nodes[tri[2]]];
        if p.iter().flatten().any(|v| !v.is_finite()) { continue; }
        let twice_area = (p[1][0]-p[0][0])*(p[2][1]-p[0][1])
            -(p[2][0]-p[0][0])*(p[1][1]-p[0][1]);
        if twice_area.abs() <= 1e-14 { continue; }

        let mut crossings = Vec::new();
        for (a,b) in [(0,1),(1,2),(2,0)] {
            let za = p[a][2]; let zb = p[b][2];
            // A strictly lower endpoint and an at-or-above endpoint.
            if (za < level) == (zb < level) { continue; }
            let t = (level - za) / (zb - za);
            if !t.is_finite() { continue; }
            crossings.push([
                p[a][0] + t*(p[b][0]-p[a][0]),
                p[a][1] + t*(p[b][1]-p[a][1]),
            ]);
        }
        if crossings.len() == 2 {
            let a = crossings[0]; let b = crossings[1];
            if a != b {
                result.push(ContourSegment { elevation: level, start: a, end: b });
            }
        }
    }
    result
}

#[cfg(test)]
mod tests {
    use super::*;
    fn tin(nodes: Vec<[f64;3]>, triangles: Vec<[usize;3]>) -> Surface {
        Surface { nodes, triangles }
    }
    #[test]
    fn inclined_triangle_interpolates_half_height() {
        let s = tin(vec![[0.,0.,0.],[2.,0.,2.],[0.,2.,0.]], vec![[0,1,2]]);
        let out = segments_at_level(&s, 1.);
        assert_eq!(out.len(), 1);
        let mut xy = [out[0].start, out[0].end];
        xy.sort_by(|a,b| a[1].total_cmp(&b[1]));
        assert!((xy[0][0]-1.).abs() < 1e-12);
        assert!((xy[0][1]-0.).abs() < 1e-12);
        assert!((xy[1][0]-1.).abs() < 1e-12);
        assert!((xy[1][1]-1.).abs() < 1e-12);
    }
    #[test]
    fn outside_range_is_empty() {
        let s = tin(vec![[0.,0.,0.],[1.,0.,1.],[0.,1.,2.]],vec![[0,1,2]]);
        assert!(segments_at_level(&s, 3.).is_empty());
    }
    #[test]
    fn flat_triangle_produces_no_spurious_segments() {
        let s = tin(vec![[0.,0.,1.],[1.,0.,1.],[0.,1.,1.]],vec![[0,1,2]]);
        assert!(segments_at_level(&s, 1.).is_empty());
    }
    #[test]
    fn degenerate_and_invalid_triangles_are_skipped() {
        let s = tin(vec![[0.,0.,0.],[1.,0.,1.],[2.,0.,2.]],vec![[0,1,2],[0,1,99]]);
        assert!(segments_at_level(&s, 1.).is_empty());
    }
    #[test]
    fn nonfinite_input_is_ignored() {
        let s = tin(vec![[0.,0.,0.],[1.,0.,f64::NAN],[0.,1.,2.]], vec![[0,1,2]]);
        assert!(segments_at_level(&s, 1.).is_empty());
    }
}
