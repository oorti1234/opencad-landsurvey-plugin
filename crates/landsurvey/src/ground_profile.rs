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
