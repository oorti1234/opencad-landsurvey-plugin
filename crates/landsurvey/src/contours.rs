//! Host-independent TIN contour sampling.
//!
//! This module reuses the isoline segments produced by the established
//! surface/datum implementation. Joining and CAD insertion are separate steps.
use crate::surface::Surface;

/// XY segment at one elevation. Coordinates use drawing units.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ContourSegment {
    pub elevation: f64,
    pub start: [f64; 2],
    pub end: [f64; 2],
}

/// Extract contour segments at one finite elevation.
///
/// Invalid triangle indices, degenerate triangles, and non-finite coordinates
/// are removed before calling the existing datum-intersection engine.
pub fn segments_at_level(surface: &Surface, elevation: f64) -> Vec<ContourSegment> {
    if !elevation.is_finite() {
        return Vec::new();
    }
    let triangles = surface.triangles.iter().copied().filter(|tri| {
        if tri.iter().any(|&index| index >= surface.nodes.len()) { return false; }
        let [a, b, c] = tri.map(|i| surface.nodes[i]);
        if [a, b, c].iter().flatten().any(|v| !v.is_finite()) { return false; }
        let twice_area = (b[0] - a[0]) * (c[1] - a[1])
            - (b[1] - a[1]) * (c[0] - a[0]);
        twice_area.is_finite() && twice_area.abs() > 1e-12
    }).collect();
    let clean = Surface { nodes: surface.nodes.clone(), triangles };
    let (_, segments) = clean.cut_fill_to_datum_detailed(elevation);
    segments.into_iter().filter_map(|[start, end]| {
        if start == end || [start, end].iter().flatten().any(|v| !v.is_finite()) {
            return None;
        }
        Some(ContourSegment { elevation, start, end })
    }).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn single_sloping_triangle_intersects_at_half_height() {
        let surface = Surface {
            nodes: vec![[0., 0., 0.], [2., 0., 2.], [0., 2., 0.]],
            triangles: vec![[0, 1, 2]],
        };
        let result = segments_at_level(&surface, 1.);
        assert_eq!(result.len(), 1);
        let mut ends = [result[0].start, result[0].end];
        ends.sort_by(|a, b| a[1].total_cmp(&b[1]));
        assert!((ends[0][0] - 1.).abs() < 1e-10);
        assert!((ends[1][0] - 1.).abs() < 1e-10);
        assert!((ends[1][1] - 1.).abs() < 1e-10);
    }

    #[test]
    fn out_of_range_and_invalid_levels_are_empty() {
        let surface = Surface {
            nodes: vec![[0., 0., 0.], [2., 0., 2.], [0., 2., 0.]],
            triangles: vec![[0, 1, 2]],
        };
        assert!(segments_at_level(&surface, 3.).is_empty());
        assert!(segments_at_level(&surface, f64::NAN).is_empty());
    }

    #[test]
    fn invalid_or_degenerate_triangles_are_skipped() {
        let surface = Surface {
            nodes: vec![[0., 0., 0.], [1., 0., 1.], [2., 0., 2.]],
            triangles: vec![[0, 1, 2], [0, 1, 99]],
        };
        assert!(segments_at_level(&surface, 1.).is_empty());
    }
}
