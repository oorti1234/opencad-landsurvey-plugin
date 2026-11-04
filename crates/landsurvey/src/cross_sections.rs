//! Exact cross-section traces through a TIN at regular axis stations.
//!
//! This is host-independent geometry. Sections use the local forward tangent
//! of the plan alignment. Offsets are **positive right**, in the direction of
//! travel; negative offsets lie to the left. At an exact axis corner, the
//! outgoing tangent is used (except at the last point). Each transverse trace
//! includes all crossed TIN edges, so terrain slope breaks are preserved.
//!
//! Do not connect samples whose `elevation` is `None` when drawing.

use crate::ground_profile::{
    trace_ground_profile, GroundProfileSample, ProfileSampleError,
};
use crate::surface::Surface;

/// Terrain sample in one cross-section. Positive offset is right of axis.
#[derive(Debug, Clone, PartialEq)]
pub struct SectionPoint {
    pub offset: f64,
    pub east: f64,
    pub north: f64,
    pub elevation: Option<f64>,
}

/// Cross-section at a centerline station, in the same linear units as TIN.
#[derive(Debug, Clone, PartialEq)]
pub struct TinCrossSection {
    pub station: f64,
    pub center: [f64; 2],
    /// Forward tangent expressed as unit vector in the XY plane.
    pub tangent: [f64; 2],
    pub points: Vec<SectionPoint>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum CrossSectionError {
    InvalidAxis,
    InvalidInterval,
    InvalidWidths,
    InvalidStartStation,
    TooManySections,
    Terrain(ProfileSampleError),
}

const MAX_SECTIONS: usize = 2048;

/// Generate sections along an XY axis at `interval` and also at its last
/// endpoint (if not an interval multiple).
///
/// - The first station is `start_station`.
/// - `left_width` and `right_width` are positive distances from centerline.
/// - A section is drawn geometrically from left through center to right.
/// - TIN-edge crossings are traced exactly, not approximated using a grid.
/// - Out-of-TIN samples retain `None` for the elevation.
/// - An axis with duplicate consecutive vertices is supported.
/// - The returned sample points are ordered by increasing offset.
pub fn trace_tin_cross_sections(
    surface: &Surface,
    axis: &[[f64; 2]],
    start_station: f64,
    interval: f64,
    left_width: f64,
    right_width: f64,
) -> Result<Vec<TinCrossSection>, CrossSectionError> {
    if !start_station.is_finite() {
        return Err(CrossSectionError::InvalidStartStation);
    }
    if !interval.is_finite() || interval <= 0.0 {
        return Err(CrossSectionError::InvalidInterval);
    }
    if !left_width.is_finite() || !right_width.is_finite()
        || left_width <= 0.0 || right_width <= 0.0
    {
        return Err(CrossSectionError::InvalidWidths);
    }
    if axis.len() < 2 || axis.iter().flatten().any(|x| !x.is_finite()) {
        return Err(CrossSectionError::InvalidAxis);
    }

    let mut segments = Vec::new();
    let mut total = 0.0;
    for pair in axis.windows(2) {
        let dx = pair[1][0] - pair[0][0];
        let dy = pair[1][1] - pair[0][1];
        let len = dx.hypot(dy);
        if !len.is_finite() {
            return Err(CrossSectionError::InvalidAxis);
        }
        if len == 0.0 {
            continue;
        }
        segments.push((total, len, pair[0], [dx / len, dy / len]));
        total += len;
        if !total.is_finite() {
            return Err(CrossSectionError::InvalidAxis);
        }
    }
    if segments.is_empty() {
        return Err(CrossSectionError::InvalidAxis);
    }

    let full = (total / interval).floor();
    if !full.is_finite() || full > MAX_SECTIONS as f64 - 2.0 {
        return Err(CrossSectionError::TooManySections);
    }
    let n = full as usize;
    let mut distances = (0..=n).map(|i| (i as f64) * interval).collect::<Vec<_>>();
    if total - *distances.last().unwrap() > 1e-9 * total.max(1.0) {
        distances.push(total);
    } else {
        *distances.last_mut().unwrap() = total;
    }
    if distances.len() > MAX_SECTIONS {
        return Err(CrossSectionError::TooManySections);
    }

    let mut sections = Vec::with_capacity(distances.len());
    for distance in distances {
        // At a vertex select the outgoing tangent, at the final endpoint
        // select the last segment's tangent.
        let &(origin_sta, segment_len, origin, tangent) = segments
            .iter()
            .find(|(s, len, _, _)| distance < s + len - 1e-10)
            .unwrap_or_else(|| segments.last().unwrap());
        let delta = (distance - origin_sta).clamp(0.0, segment_len);
        let center = [
            origin[0] + delta * tangent[0],
            origin[1] + delta * tangent[1],
        ];

        // Right normal for a forward tangent (dx,dy) is (dy,-dx).
        let right = [tangent[1], -tangent[0]];
        let left_pt = [
            center[0] - left_width * right[0],
            center[1] - left_width * right[1],
        ];
        let right_pt = [
            center[0] + right_width * right[0],
            center[1] + right_width * right[1],
        ];
        if left_pt.iter().chain(right_pt.iter()).any(|x| !x.is_finite()) {
            return Err(CrossSectionError::InvalidWidths);
        }

        // Including the centerline as an axis vertex guarantees an offset=0
        // sample even when no TIN edge happens to cross the centerline.
        let profile: Vec<GroundProfileSample> = trace_ground_profile(
            surface,
            &[left_pt, center, right_pt],
            -left_width,
        )
        .map_err(CrossSectionError::Terrain)?;
        let points = profile.into_iter().map(|p| SectionPoint {
            offset: p.station,
            east: p.east,
            north: p.north,
            elevation: p.elevation,
        }).collect();
        let station = start_station + distance;
        if !station.is_finite() {
            return Err(CrossSectionError::InvalidStartStation);
        }
        sections.push(TinCrossSection { station, center, tangent, points });
    }
    Ok(sections)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane() -> Surface {
        // TIN square: Z=100+E+2*N.
        Surface {
            nodes: vec![
                [0., 0., 100.], [10., 0., 110.],
                [0., 10., 120.], [10., 10., 130.],
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        }
    }

    #[test]
    fn east_axis_offsets_positive_right_and_elevations_exact() {
        let ss = trace_tin_cross_sections(
            &plane(), &[[2., 5.], [8., 5.]],
            100., 3., 2., 2.,
        ).unwrap();
        assert_eq!(ss.len(), 3);
        assert_eq!(ss.iter().map(|s| s.station).collect::<Vec<_>>(),
                   vec![100., 103., 106.]);
        let s = &ss[1];
        assert_eq!(s.center, [5., 5.]);
        assert_eq!(s.tangent, [1., 0.]);
        let left = s.points.first().unwrap();
        let right = s.points.last().unwrap();
        assert!((left.offset + 2.).abs() < 1e-9);
        assert!((right.offset - 2.).abs() < 1e-9);
        assert!((left.north - 7.).abs() < 1e-9);
        assert!((right.north - 3.).abs() < 1e-9);
        assert!((left.elevation.unwrap() - 119.).abs() < 1e-9);
        assert!((right.elevation.unwrap() - 111.).abs() < 1e-9);
        assert!(s.points.iter().any(|p| p.offset.abs() < 1e-9
            && (p.elevation.unwrap() - 115.).abs() < 1e-9));
    }

    #[test]
    fn last_station_is_included_even_when_interval_is_not_a_divisor() {
        let ss = trace_tin_cross_sections(
            &plane(), &[[1., 5.], [9., 5.]],
            0., 3., 1., 1.,
        ).unwrap();
        assert_eq!(ss.iter().map(|s| s.station).collect::<Vec<_>>(),
                   vec![0., 3., 6., 8.]);
    }

    #[test]
    fn section_at_corner_uses_outgoing_direction() {
        let ss = trace_tin_cross_sections(
            &plane(), &[[2., 2.], [6., 2.], [6., 6.]],
            0., 4., 1., 1.,
        ).unwrap();
        assert_eq!(ss.len(), 3);
        assert_eq!(ss[1].center, [6., 2.]);
        assert_eq!(ss[1].tangent, [0., 1.]);
        // Forward north means right/east, hence offset+1 is east of axis.
        assert!((ss[1].points.last().unwrap().east - 7.).abs() < 1e-9);
    }

    #[test]
    fn outside_tin_does_not_extrapolate() {
        let ss = trace_tin_cross_sections(
            &plane(), &[[5., 2.], [5., 8.]],
            0., 6., 10., 10.,
        ).unwrap();
        assert!(ss[0].points.first().unwrap().elevation.is_none());
        assert!(ss[0].points.last().unwrap().elevation.is_none());
        assert!(ss[0].points.iter().any(|p| p.elevation.is_some()));
    }

    #[test]
    fn invalid_parameters_and_empty_axis_are_rejected() {
        assert_eq!(trace_tin_cross_sections(&plane(), &[[0.,0.],[1.,0.]], 0., 0., 2., 2.),
                   Err(CrossSectionError::InvalidInterval));
        assert_eq!(trace_tin_cross_sections(&plane(), &[[0.,0.],[1.,0.]], 0., 1., -2., 2.),
                   Err(CrossSectionError::InvalidWidths));
        assert_eq!(trace_tin_cross_sections(&plane(), &[[0.,0.],[0.,0.]], 0., 1., 2., 2.),
                   Err(CrossSectionError::InvalidAxis));
        assert_eq!(trace_tin_cross_sections(&plane(), &[[0.,0.],[100000.,0.]], 0., 0.001, 2., 2.),
                   Err(CrossSectionError::TooManySections));
    }
}
