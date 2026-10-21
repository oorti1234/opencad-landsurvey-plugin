//! Longitudinal terrain sampling along an XY alignment.
//! Stations are horizontal distances from the first vertex. Sampling uses
//! the existing indexed TIN interpolation, without changing surface geometry.
use crate::surface::Surface;
use crate::tin_spatial_index::TinGridIndex;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ProfileSample {
    pub station: f64,
    pub x: f64,
    pub y: f64,
    /// None means the sample falls outside the triangulated surface.
    pub elevation: Option<f64>,
}

#[derive(Debug, Clone, PartialEq)]
pub struct LongitudinalProfile {
    pub length: f64,
    pub samples: Vec<ProfileSample>,
    pub outside_count: usize,
}

pub fn sample_longitudinal(
    surface: &Surface,
    alignment: &[[f64; 2]],
    spacing: f64,
) -> Result<LongitudinalProfile, String> {
    if !spacing.is_finite() || spacing <= 0.0 {
        return Err("Profile spacing must be finite and positive".into());
    }
    if alignment.len() < 2 || alignment.iter().flatten().any(|v| !v.is_finite()) {
        return Err("Profile requires at least two finite XY vertices".into());
    }
    let mut cumulative = Vec::with_capacity(alignment.len());
    cumulative.push(0.0);
    for segment in alignment.windows(2) {
        let d = (segment[1][0] - segment[0][0]).hypot(segment[1][1] - segment[0][1]);
        if !d.is_finite() || d <= 0.0 {
            return Err("Profile alignment has a zero-length or invalid segment".into());
        }
        let next = cumulative.last().copied().unwrap_or(0.0) + d;
        if !next.is_finite() {
            return Err("Profile alignment length overflow".into());
        }
        cumulative.push(next);
    }
    let length = *cumulative.last().unwrap();
    let interior = (length / spacing).ceil();
    if !interior.is_finite() || interior > 200_000.0 {
        return Err("Profile has too many sample stations (limit 200000)".into());
    }
    let index = TinGridIndex::build(surface)?;
    let steps = (length / spacing).floor() as usize;
    let mut stations: Vec<f64> = (0..=steps).map(|k| k as f64 * spacing).collect();
    if let Some(last) = stations.last_mut() {
        if (length - *last).abs() <= length.max(1.0) * 1e-12 {
            *last = length;
        } else {
            stations.push(length);
        }
    }
    let mut samples = Vec::with_capacity(stations.len());
    let mut seg = 0usize;
    let mut outside_count = 0;
    for station in stations {
        while seg + 1 < cumulative.len() - 1 && station > cumulative[seg + 1] {
            seg += 1;
        }
        let t = ((station - cumulative[seg]) / (cumulative[seg + 1] - cumulative[seg])).clamp(0.0, 1.0);
        let x = alignment[seg][0] + t * (alignment[seg + 1][0] - alignment[seg][0]);
        let y = alignment[seg][1] + t * (alignment[seg + 1][1] - alignment[seg][1]);
        let elevation = index.interpolate_z(surface, x, y);
        if elevation.is_none() { outside_count += 1; }
        samples.push(ProfileSample { station, x, y, elevation });
    }
    Ok(LongitudinalProfile { length, samples, outside_count })
}

#[cfg(test)]
mod tests {
    use super::*;
    fn planar() -> Surface {
        Surface { nodes: vec![[0.0,0.0,100.0],[10.0,0.0,110.0],
                              [0.0,10.0,120.0],[10.0,10.0,130.0]],
                  triangles: vec![[0,1,2],[1,3,2]] }
    }
    #[test]
    fn samples_plane_at_fixed_stations_and_endpoint() {
        let p = sample_longitudinal(&planar(), &[[0.0,0.0],[10.0,0.0]], 3.0).unwrap();
        assert_eq!(p.samples.len(),5);
        assert_eq!(p.outside_count,0);
        for (i,s) in p.samples.iter().enumerate() {
            let expected = [0.0,3.0,6.0,9.0,10.0][i];
            assert!((s.station-expected).abs()<1e-9);
            assert!((s.elevation.unwrap()-(100.0+expected)).abs()<1e-9);
        }
    }
    #[test]
    fn polyline_stations_follow_cumulative_length() {
        let p=sample_longitudinal(&planar(), &[[0.0,0.0],[3.0,0.0],[3.0,4.0]], 2.0).unwrap();
        assert_eq!(p.samples.iter().map(|s|s.station).collect::<Vec<_>>(), vec![0.0,2.0,4.0,6.0,7.0]);
        assert!((p.samples[2].x-3.0).abs()<1e-9);
        assert!((p.samples[2].y-1.0).abs()<1e-9);
        assert!((p.samples[2].elevation.unwrap()-105.0).abs()<1e-9);
    }
    #[test]
    fn marks_outside_without_inventing_elevation() {
        let p=sample_longitudinal(&planar(), &[[-2.0,0.0],[2.0,0.0]], 2.0).unwrap();
        assert_eq!(p.outside_count,1);
        assert_eq!(p.samples[0].elevation,None);
        assert!((p.samples[2].elevation.unwrap()-102.0).abs()<1e-9);
    }
    #[test]
    fn rejects_bad_spacing_and_duplicate_vertices() {
        assert!(sample_longitudinal(&planar(), &[[0.0,0.0],[2.0,0.0]], 0.0).is_err());
        assert!(sample_longitudinal(&planar(), &[[0.0,0.0],[0.0,0.0]], 1.0).is_err());
        assert!(sample_longitudinal(&planar(), &[[0.0,0.0],[2.0,0.0]], 1e-9).is_err());
    }
}


/// Design-grade elevation and signed difference at each terrain profile station.
/// Positive difference means planned fill; negative means planned cut.
/// These are vertical differences, NOT earthwork volumes.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GradeSample {
    pub station: f64,
    pub terrain_z: Option<f64>,
    pub design_z: f64,
    pub design_minus_terrain: Option<f64>,
}

/// Evaluate a constant longitudinal slope in percent from the first station.
/// The design equation is Z(s) = initial_z + s * slope_percent / 100.
pub fn constant_grade(
    profile: &LongitudinalProfile,
    initial_z: f64,
    slope_percent: f64,
) -> Result<Vec<GradeSample>, String> {
    if !initial_z.is_finite() || !slope_percent.is_finite() {
        return Err("Grade elevation and slope must be finite".into());
    }
    if !profile.length.is_finite() || profile.length < 0.0 {
        return Err("Invalid profile length".into());
    }
    let mut last_station = -1.0;
    let mut result = Vec::with_capacity(profile.samples.len());
    for sample in &profile.samples {
        if !sample.station.is_finite() || sample.station < 0.0
            || sample.station > profile.length + 1e-8 || sample.station <= last_station {
            return Err("Invalid or unordered profile station".into());
        }
        last_station = sample.station;
        let design_z = initial_z + sample.station * (slope_percent / 100.0);
        if !design_z.is_finite() {
            return Err("Grade elevation overflow".into());
        }
        let difference = match sample.elevation {
            Some(z) if z.is_finite() => Some(design_z - z),
            Some(_) => return Err("Invalid terrain elevation".into()),
            None => None,
        };
        result.push(GradeSample {
            station: sample.station,
            terrain_z: sample.elevation,
            design_z,
            design_minus_terrain: difference,
        });
    }
    Ok(result)
}

#[cfg(test)]
mod grade_tests {
    use super::*;

    #[test]
    fn rise_and_signed_cut_fill_are_consistent() {
        let ground = sample_longitudinal(
            &Surface { nodes: vec![[0.,0.,100.],[10.,0.,110.],[0.,10.,100.]],
                triangles: vec![[0,1,2]] },
            &[[0.,0.],[10.,0.]], 5.).unwrap();
        let values = constant_grade(&ground, 102., 20.).unwrap();
        assert_eq!(values.len(), 3);
        assert!((values[0].design_minus_terrain.unwrap()-2.).abs()<1e-9);
        assert!((values[1].design_z-103.).abs()<1e-9);
        assert!((values[1].design_minus_terrain.unwrap()+2.).abs()<1e-9);
        assert!((values[2].design_minus_terrain.unwrap()+6.).abs()<1e-9);
    }

    #[test]
    fn falling_grade_and_missing_terrain() {
        let p = LongitudinalProfile { length:10., outside_count:1,
            samples:vec![
                ProfileSample{station:0.,x:0.,y:0.,elevation:Some(104.)},
                ProfileSample{station:5.,x:5.,y:0.,elevation:None},
                ProfileSample{station:10.,x:10.,y:0.,elevation:Some(102.)}
            ]};
        let grade=constant_grade(&p, 103., -10.).unwrap();
        assert_eq!(grade[0].design_minus_terrain, Some(-1.));
        assert_eq!(grade[1].design_minus_terrain,None);
        assert_eq!(grade[2].design_minus_terrain,Some(0.));
    }

    #[test]
    fn rejects_invalid_input() {
        let p=LongitudinalProfile { length:5.,outside_count:0,samples:vec![
            ProfileSample{station:0.,x:0.,y:0.,elevation:Some(100.)}
        ]};
        assert!(constant_grade(&p,f64::NAN,1.).is_err());
        assert!(constant_grade(&p,100.,f64::INFINITY).is_err());
    }
}
