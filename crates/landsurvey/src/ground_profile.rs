//! Longitudinal ground-profile sampling over a TIN (host-independent).
//!
//! Input is a plan-view polyline in [Easting, Northing] order. Samples are
//! spaced no farther than the requested step *within each axis segment*;
//! original corners are always retained. This is a sampled profile, not an
//! exact TIN-face intersection profile: short grade breaks may lie between
//! sample stations. A missing elevation means that the station is outside
//! the TIN; callers must not draw a continuous line through such gaps.

use crate::surface::Surface;

/// One alignment station and the TIN elevation (if within the TIN footprint).
#[derive(Debug, Clone, PartialEq)]
pub struct GroundProfileSample {
    pub station: f64,
    pub east: f64,
    pub north: f64,
    pub elevation: Option<f64>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ProfileSampleError {
    InvalidStep,
    InvalidStartStation,
    InvalidAxis,
    InvalidSurface,
    TooManySamples,
}

const MAX_SAMPLES: usize = 100_000;

/// Sample a longitudinal profile along a plan-view polyline.
///
/// - `step` is a positive maximum sample distance (same units as the TIN).
/// - `start_station` is the first station (e.g. 1000.0 for 1+000).
/// - The axis must have at least two distinct points; coincident consecutive
///   axis points are ignored without creating duplicate stations.
/// - Each axis corner is sampled exactly once.
/// - At a station outside the TIN, `elevation` is `None`.
/// - This API does not smooth terrain, extrapolate outside the TIN, or return
///   CAD entities. Rendering and XDATA belong in the plugin layer.
pub fn sample_ground_profile(
    surface: &Surface,
    axis: &[[f64; 2]],
    start_station: f64,
    step: f64,
) -> Result<Vec<GroundProfileSample>, ProfileSampleError> {
    if !step.is_finite() || step <= 0.0 {
        return Err(ProfileSampleError::InvalidStep);
    }
    if !start_station.is_finite() {
        return Err(ProfileSampleError::InvalidStartStation);
    }
    if axis.len() < 2 || axis.iter().flatten().any(|v| !v.is_finite()) {
        return Err(ProfileSampleError::InvalidAxis);
    }
    if surface.nodes.is_empty()
        || surface.triangles.is_empty()
        || surface.nodes.iter().flatten().any(|v| !v.is_finite())
        || surface.triangles.iter().flatten().any(|&i| i >= surface.nodes.len())
    {
        return Err(ProfileSampleError::InvalidSurface);
    }

    // First calculate the required number of samples; avoid allocating an
    // unbounded vector on malformed or extremely long axes.
    let mut required = 1usize;
    let mut valid_segments = 0usize;
    for w in axis.windows(2) {
        let dx = w[1][0] - w[0][0];
        let dy = w[1][1] - w[0][1];
        let length = dx.hypot(dy);
        if !length.is_finite() {
            return Err(ProfileSampleError::InvalidAxis);
        }
        if length == 0.0 {
            continue;
        }
        let divisions = (length / step).ceil();
        if !divisions.is_finite() || divisions > MAX_SAMPLES as f64 {
            return Err(ProfileSampleError::TooManySamples);
        }
        let n = (divisions as usize).max(1);
        required = required.checked_add(n).ok_or(ProfileSampleError::TooManySamples)?;
        if required > MAX_SAMPLES {
            return Err(ProfileSampleError::TooManySamples);
        }
        valid_segments += 1;
    }
    if valid_segments == 0 {
        return Err(ProfileSampleError::InvalidAxis);
    }

    let mut samples = Vec::with_capacity(required);
    let mut traveled = 0.0f64;
    let mut first = true;
    for w in axis.windows(2) {
        let dx = w[1][0] - w[0][0];
        let dy = w[1][1] - w[0][1];
        let length = dx.hypot(dy);
        if length == 0.0 {
            continue;
        }
        let count = ((length / step).ceil() as usize).max(1);
        if first {
            let p = w[0];
            samples.push(GroundProfileSample {
                station: start_station + traveled,
                east: p[0],
                north: p[1],
                elevation: surface.interpolate_z(p[0], p[1]),
            });
            first = false;
        }
        for j in 1..=count {
            let t = j as f64 / count as f64;
            let east = w[0][0] + t * dx;
            let north = w[0][1] + t * dy;
            let station = start_station + traveled + t * length;
            if !station.is_finite() {
                return Err(ProfileSampleError::InvalidAxis);
            }
            samples.push(GroundProfileSample {
                station,
                east,
                north,
                elevation: surface.interpolate_z(east, north),
            });
        }
        traveled += length;
        if !traveled.is_finite() {
            return Err(ProfileSampleError::InvalidAxis);
        }
    }
    Ok(samples)
}

/// Trace a longitudinal ground profile at every intersection with a TIN edge.
///
/// Unlike fixed-step sampling, this preserves the exact slope breaks of the
/// piecewise-planar TIN along each straight segment of the supplied axis.
/// Interior TIN edges and boundary edges are considered. Where the axis runs
/// outside the surface footprint, a sample with `elevation: None` separates
/// drawable spans. Collinear boundary overlaps are handled by including the
/// projected edge endpoints. The output includes the axis start/end and
/// corners. No CAD entities are created.
///
/// This is intended for valid, non-overlapping TIN triangles. Geometry at
/// ambiguous coincident or overlapping triangles is resolved by the existing
/// Surface::interpolate_z() policy.
pub fn trace_ground_profile(
    surface: &Surface,
    axis: &[[f64; 2]],
    start_station: f64,
) -> Result<Vec<GroundProfileSample>, ProfileSampleError> {
    if !start_station.is_finite() {
        return Err(ProfileSampleError::InvalidStartStation);
    }
    if axis.len() < 2 || axis.iter().flatten().any(|v| !v.is_finite()) {
        return Err(ProfileSampleError::InvalidAxis);
    }
    if surface.nodes.is_empty() || surface.triangles.is_empty()
        || surface.nodes.iter().flatten().any(|v| !v.is_finite())
        || surface.triangles.iter().flatten().any(|&i| i >= surface.nodes.len())
    {
        return Err(ProfileSampleError::InvalidSurface);
    }

    let cross = |a: [f64; 2], b: [f64; 2]| a[0] * b[1] - a[1] * b[0];
    let mut result = Vec::new();
    let mut traveled = 0.0;
    let mut any_segment = false;
    for segment in axis.windows(2) {
        let p = segment[0];
        let r = [segment[1][0] - p[0], segment[1][1] - p[1]];
        let len = r[0].hypot(r[1]);
        if !len.is_finite() {
            return Err(ProfileSampleError::InvalidAxis);
        }
        if len == 0.0 {
            continue;
        }
        any_segment = true;
        // Parameter space is local to the axis segment, [0, 1].
        let mut ts = vec![0.0f64, 1.0];
        for &tri in &surface.triangles {
            for (ia, ib) in [(tri[0], tri[1]), (tri[1], tri[2]), (tri[2], tri[0])] {
                let a = surface.nodes[ia];
                let b = surface.nodes[ib];
                let q = [a[0], a[1]];
                let edge = [b[0] - q[0], b[1] - q[1]];
                let qp = [q[0] - p[0], q[1] - p[1]];
                let den = cross(r, edge);
                // Scale the cross-product tolerance to the local dimensions.
                let tol = 1e-12 * len * edge[0].hypot(edge[1]);
                if den.abs() > tol {
                    let t = cross(qp, edge) / den;
                    let u = cross(qp, r) / den;
                    if (-1e-10..=1.0 + 1e-10).contains(&t)
                        && (-1e-10..=1.0 + 1e-10).contains(&u)
                    {
                        ts.push(t.clamp(0.0, 1.0));
                    }
                } else if cross(qp, r).abs() <= 1e-12 * len * qp[0].hypot(qp[1]).max(1.0) {
                    // Collinear axis/edge: include both edge endpoints.
                    let rr = r[0] * r[0] + r[1] * r[1];
                    let t0 = (qp[0] * r[0] + qp[1] * r[1]) / rr;
                    let qb = [b[0] - p[0], b[1] - p[1]];
                    let t1 = (qb[0] * r[0] + qb[1] * r[1]) / rr;
                    for t in [t0, t1] {
                        if (-1e-10..=1.0 + 1e-10).contains(&t) {
                            ts.push(t.clamp(0.0, 1.0));
                        }
                    }
                }
            }
        }
        ts.sort_by(|a, b| a.total_cmp(b));
        ts.dedup_by(|a, b| (*a - *b).abs() <= 1e-10);
        // Count both station breakpoints and possible gap markers.
        if result.len().saturating_add(ts.len().saturating_mul(2)) > MAX_SAMPLES {
            return Err(ProfileSampleError::TooManySamples);
        }
        let point_at = |t: f64| [p[0] + t * r[0], p[1] + t * r[1]];
        for (i, &t) in ts.iter().enumerate() {
            if !result.is_empty() && i == 0 {
                continue; // axis corner already in the output
            }
            if i > 0 {
                let prev = ts[i - 1];
                let mid = (prev + t) * 0.5;
                let m = point_at(mid);
                if surface.interpolate_z(m[0], m[1]).is_none() {
                    result.push(GroundProfileSample {
                        station: start_station + traveled + mid * len,
                        east: m[0],
                        north: m[1],
                        elevation: None,
                    });
                }
            }
            let xy = point_at(t);
            let station = start_station + traveled + t * len;
            if !station.is_finite() {
                return Err(ProfileSampleError::InvalidAxis);
            }
            result.push(GroundProfileSample {
                station,
                east: xy[0],
                north: xy[1],
                elevation: surface.interpolate_z(xy[0], xy[1]),
            });
        }
        traveled += len;
        if !traveled.is_finite() {
            return Err(ProfileSampleError::InvalidAxis);
        }
    }
    if !any_segment {
        return Err(ProfileSampleError::InvalidAxis);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn plane() -> Surface {
        // Square with Z = 100 + E + 2*N.
        Surface {
            nodes: vec![
                [0., 0., 100.], [10., 0., 110.],
                [0., 10., 120.], [10., 10., 130.],
            ],
            triangles: vec![[0, 1, 2], [1, 3, 2]],
        }
    }

    #[test]
    fn exact_trace_keeps_tin_edge_grade_break() {
        // Two adjacent facets with a sharp change in longitudinal slope at E=5.
        let tin = Surface {
            nodes: vec![
                [0., 0., 100.], [5., 0., 105.], [10., 0., 100.],
                [0., 10., 100.], [5., 10., 105.], [10., 10., 100.],
            ],
            triangles: vec![[0, 1, 3], [1, 4, 3], [1, 2, 4], [2, 5, 4]],
        };
        let p = trace_ground_profile(&tin, &[[0., 5.], [10., 5.]], 1000.).unwrap();
        assert!(p.iter().any(|x|
            (x.station - 1005.).abs() < 1e-9
            && (x.elevation.unwrap() - 105.).abs() < 1e-9
        ), "TIN break at station 1+005 absent: {p:?}");
        assert_eq!(p.first().unwrap().elevation, Some(100.));
        assert_eq!(p.last().unwrap().elevation, Some(100.));
    }

    #[test]
    fn exact_trace_marks_outside_gaps() {
        let p = trace_ground_profile(
            &plane(), &[[-5., 5.], [15., 5.]], 0.,
        ).unwrap();
        assert_eq!(p.first().unwrap().elevation, None);
        assert_eq!(p.last().unwrap().elevation, None);
        assert!(p.iter().any(|x| x.elevation.is_some()));
        assert!(p.iter().any(|x| x.elevation.is_none() && x.station > 0. && x.station < 5.));
        assert!(p.iter().any(|x| x.elevation.is_none() && x.station > 15. && x.station < 20.));
    }

    #[test]
    fn exact_trace_keeps_polyline_corners_and_stations() {
        let p = trace_ground_profile(
            &plane(), &[[1., 1.], [4., 1.], [4., 4.]], 50.,
        ).unwrap();
        assert!(p.iter().any(|x| (x.station - 53.).abs() < 1e-10
            && (x.east - 4.).abs() < 1e-10 && (x.north - 1.).abs() < 1e-10));
        assert!((p.last().unwrap().station - 56.).abs() < 1e-10);
    }

    #[test]
    fn exact_trace_rejects_bad_axis_and_tin() {
        assert_eq!(trace_ground_profile(&plane(), &[[1., 1.]], 0.),
                   Err(ProfileSampleError::InvalidAxis));
        assert_eq!(trace_ground_profile(&plane(), &[[1., 1.], [1., 1.]], 0.),
                   Err(ProfileSampleError::InvalidAxis));
        assert_eq!(trace_ground_profile(&plane(), &[[1., 1.], [2., 1.]], f64::NAN),
                   Err(ProfileSampleError::InvalidStartStation));
    }

    #[test]
    fn known_plane_profile_and_stationing() {
        let p = sample_ground_profile(
            &plane(), &[[1., 2.], [9., 2.]], 1000., 2.,
        ).unwrap();
        assert_eq!(p.len(), 5);
        for (i, s) in p.iter().enumerate() {
            let east = 1. + 2. * i as f64;
            assert!((s.station - (1000. + 2. * i as f64)).abs() < 1e-10);
            assert!((s.east - east).abs() < 1e-10);
            assert!((s.north - 2.).abs() < 1e-10);
            assert!((s.elevation.unwrap() - (104. + east)).abs() < 1e-9);
        }
    }

    #[test]
    fn corners_are_preserved_without_duplicate_stations() {
        let p = sample_ground_profile(
            &plane(),
            &[[0., 0.], [4., 0.], [4., 0.], [4., 3.]],
            0., 2.,
        ).unwrap();
        let stations: Vec<f64> = p.iter().map(|v| v.station).collect();
        assert_eq!(stations, vec![0., 2., 4., 5.5, 7.]);
        assert_eq!((p[2].east, p[2].north), (4., 0.));
        assert_eq!((p[4].east, p[4].north), (4., 3.));
    }

    #[test]
    fn outside_tin_is_explicitly_missing_not_extrapolated() {
        let p = sample_ground_profile(
            &plane(), &[[-5., 5.], [15., 5.]], 0., 5.,
        ).unwrap();
        assert_eq!(p.len(), 5);
        assert_eq!(p.iter().map(|x| x.elevation.is_some()).collect::<Vec<_>>(),
                   vec![false, true, true, true, false]);
    }

    #[test]
    fn invalid_inputs_are_rejected() {
        assert_eq!(sample_ground_profile(&plane(), &[[0., 0.], [2., 0.]], 0., 0.),
                   Err(ProfileSampleError::InvalidStep));
        assert_eq!(sample_ground_profile(&plane(), &[[0., 0.], [2., 0.]], f64::NAN, 1.),
                   Err(ProfileSampleError::InvalidStartStation));
        assert_eq!(sample_ground_profile(&plane(), &[[0., 0.]], 0., 1.),
                   Err(ProfileSampleError::InvalidAxis));
        assert_eq!(sample_ground_profile(&plane(), &[[1., 1.], [1., 1.]], 0., 1.),
                   Err(ProfileSampleError::InvalidAxis));
        assert_eq!(sample_ground_profile(&plane(), &[[0., 0.], [1_000_000., 0.]], 0., 0.001),
                   Err(ProfileSampleError::TooManySamples));
    }

    #[test]
    fn degenerate_tin_and_nonfinite_axis_are_rejected() {
        let empty = Surface { nodes: vec![], triangles: vec![] };
        assert_eq!(sample_ground_profile(&empty, &[[0., 0.], [1., 0.]], 0., 1.),
                   Err(ProfileSampleError::InvalidSurface));
        assert_eq!(sample_ground_profile(&plane(), &[[0., 0.], [f64::INFINITY, 0.]], 0., 1.),
                   Err(ProfileSampleError::InvalidAxis));
    }
}
