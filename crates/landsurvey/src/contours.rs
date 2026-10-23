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

    /// Analytic reference independent of Civil 3D:
    /// on the plane z = x a level z = 1 is the straight line
    /// x = 1, y in [0, 2], with exact length 2 units.
    #[test]
    fn planar_slope_two_triangles_match_analytic_length() {
        let surface = Surface {
            nodes: vec![
                [0., 0., 0.], [2., 0., 2.],
                [0., 2., 0.], [2., 2., 2.],
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        };
        let segments = segments_at_level(&surface, 1.);
        assert_eq!(segments.len(), 2, "one segment per planar facet");
        let mut length = 0.0;
        for segment in segments {
            for end in [segment.start, segment.end] {
                assert!((end[0] - 1.).abs() < 1e-10, "x must equal contour elevation");
                assert!(end[1] >= -1e-10 && end[1] <= 2. + 1e-10);
            }
            length += (segment.end[0] - segment.start[0])
                .hypot(segment.end[1] - segment.start[1]);
        }
        assert!((length - 2.).abs() < 1e-10,
            "analytic reference: total contour length must be 2");
    }

    /// Analytic end-to-end reference across a shared TIN edge:
    /// the contour x=1 joins two triangle segments into one open polyline,
    /// with three vertices and total length exactly 2.
    #[test]
    fn planar_slope_stitches_into_one_analytic_contour() {
        use crate::contour_chains::stitch_segments;
        let surface = Surface {
            nodes: vec![
                [0., 0., 0.], [2., 0., 2.],
                [0., 2., 0.], [2., 2., 2.],
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        };
        let segments = segments_at_level(&surface, 1.);
        let chains = stitch_segments(&segments, 1e-8).expect("planar slope is unambiguous");
        assert_eq!(chains.len(), 1);
        assert!(!chains[0].closed);
        assert_eq!(chains[0].points.len(), 3);
        assert!(chains[0].points.iter().all(|p| (p[0] - 1.).abs() < 1e-10));
        let length: f64 = chains[0].points.windows(2)
            .map(|pair| (pair[1][0] - pair[0][0]).hypot(pair[1][1] - pair[0][1]))
            .sum();
        assert!((length - 2.).abs() < 1e-10);
    }

    /// Translation in Easting/Northing and a constant vertical offset must
    /// translate the contour without changing the total geometry length.
    #[test]
    fn planar_slope_translation_preserves_contour_length() {
        let surface = Surface {
            nodes: vec![
                [500., 1000., 20.], [502., 1000., 22.],
                [500., 1002., 20.], [502., 1002., 22.],
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        };
        let segments = segments_at_level(&surface, 21.);
        assert_eq!(segments.len(), 2);
        let length: f64 = segments.iter().map(|s| {
            (s.end[0] - s.start[0]).hypot(s.end[1] - s.start[1])
        }).sum();
        assert!((length - 2.).abs() < 1e-10);
        assert!(segments.iter().all(|s|
            [s.start, s.end].iter().all(|p| (p[0] - 501.).abs() < 1e-10)));
    }

    /// A true TIN saddle: four facets meet at a vertex exactly at level zero,
    /// with alternating positive and negative surrounding vertices.
    /// Do not guess how its four branches should connect.
    #[test]
    fn exact_vertex_saddle_is_not_silently_cross_connected() {
        use crate::contour_chains::{stitch_segments, ChainError};
        let surface = Surface {
            nodes: vec![
                [0., 0., 0.],
                [1., 0., 1.],
                [0., 1., -1.],
                [-1., 0., 1.],
                [0., -1., -1.],
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
        };
        let segments = segments_at_level(&surface, 0.);
        assert!(!segments.is_empty(), "saddle must produce nontrivial contour arms");
        assert!(matches!(
            stitch_segments(&segments, 1e-9),
            Err(ChainError::Branch { degree: 4, .. })
        ), "ambiguous four-arm saddle must be rejected, not cross-connected");
    }

    /// Moving the contour off an exact saddle must give ordinary,
    /// unambiguous paths rather than a spurious four-way intersection.
    /// The piecewise-linear surface has z=+1 on E/W and z=-1 on N/S.
    /// Above the saddle, isolines occur in the east and west wedges;
    /// below it, they occur in the north and south wedges.
    #[test]
    fn saddle_offset_levels_form_two_separate_open_chains() {
        use crate::contour_chains::stitch_segments;
        let surface = Surface {
            nodes: vec![
                [0., 0., 0.], [1., 0., 1.], [0., 1., -1.],
                [-1., 0., 1.], [0., -1., -1.],
            ],
            triangles: vec![[0, 1, 2], [0, 2, 3], [0, 3, 4], [0, 4, 1]],
        };
        for (level, on_east_west) in [(0.25, true), (-0.25, false)] {
            let segments = segments_at_level(&surface, level);
            let chains = stitch_segments(&segments, 1e-9)
                .expect("off-saddle isolines must be unambiguous");
            assert_eq!(chains.len(), 2, "expected two disconnected branches");
            assert!(chains.iter().all(|c| !c.closed && c.points.len() >= 2));
            for chain in &chains {
                let midpoint = chain.points[chain.points.len() / 2];
                if on_east_west {
                    assert!(midpoint[0].abs() > 1e-8,
                        "positive-level branches belong to east/west wedges");
                } else {
                    assert!(midpoint[1].abs() > 1e-8,
                        "negative-level branches belong to north/south wedges");
                }
            }
        }
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
