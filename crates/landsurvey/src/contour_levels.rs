//! Experimental batch contours for an existing TIN surface.
//! Geometry stays host-independent; DXF/ribbon integration comes later.
use crate::contour_chains::{stitch_segments_split_junctions, ContourChain};
use crate::contours::segments_at_level;
use crate::surface::Surface;

#[derive(Debug, Clone, PartialEq)]
pub struct ClassifiedContour {
    pub chain: ContourChain,
    pub major: bool,
}

/// Compute contours for a TIN. All elevations and distances use drawing units.
///
/// `minor_interval` must be positive; `major_every` is the number of minor
/// intervals between major contours. `base_elevation` anchors the levels.
/// `tolerance` is the XY joining tolerance in drawing units.
/// Fail closed on excessive level counts or malformed settings.
pub fn generate_contours(
    surface: &Surface,
    minor_interval: f64,
    major_every: usize,
    base_elevation: f64,
    tolerance: f64,
) -> Result<Vec<ClassifiedContour>, String> {
    if !minor_interval.is_finite() || minor_interval <= 0.0
        || !base_elevation.is_finite()
        || !tolerance.is_finite() || tolerance <= 0.0
        || major_every == 0 {
        return Err("invalid contour settings".into());
    }
    if surface.nodes.is_empty() || surface.triangles.is_empty() {
        return Ok(Vec::new());
    }
    if surface.nodes.iter().flatten().any(|v| !v.is_finite()) {
        return Err("surface contains non-finite coordinates".into());
    }
    if surface.triangles.iter().flatten().any(|i| *i >= surface.nodes.len()) {
        return Err("TIN contains out-of-range vertex indices".into());
    }
    let lo = surface.nodes.iter().map(|p|p[2]).fold(f64::INFINITY, f64::min);
    let hi = surface.nodes.iter().map(|p|p[2]).fold(f64::NEG_INFINITY, f64::max);
    let first = ((lo - base_elevation)/minor_interval).ceil();
    let last = ((hi - base_elevation)/minor_interval).floor();
    if !first.is_finite() || !last.is_finite()
        || first < i64::MIN as f64 || last >= i64::MAX as f64 {
        return Err("contour index outside supported range".into());
    }
    if first > last { return Ok(Vec::new()); }
    if last - first > 10000.0 {
        return Err("too many contour levels (maximum 10001)".into());
    }
    let mut result = Vec::new();
    for n in (first as i64)..=(last as i64) {
        let level = base_elevation + (n as f64)*minor_interval;
        if !level.is_finite() { return Err("non-finite contour level".into()); }
        let segments = segments_at_level(surface, level);
        let chains = stitch_segments_split_junctions(&segments, tolerance)?;
        for chain in chains {
            result.push(ClassifiedContour {
                chain, major: n.rem_euclid(major_every as i64) == 0,
            });
        }
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn ramp() -> Surface {
        Surface {
            nodes: vec![[0.,0.,0.],[10.,0.,5.],[0.,10.,0.]],
            triangles: vec![[0,1,2]],
        }
    }
    #[test]
    fn generates_half_meter_levels_and_2_5m_majors() {
        let out=generate_contours(&ramp(),0.5,5,0.,1e-6).unwrap();
        assert_eq!(out.len(),9); // 0.5 through 4.5; endpoints do not yield segments
        assert_eq!(out.iter().filter(|c|c.major).count(),1);
        assert!((out.iter().find(|c|c.major).unwrap().chain.elevation-2.5).abs()<1e-12);
    }
    #[test]
    fn level_vertex_with_four_branches_does_not_abort_generation() {
        let s=Surface{nodes:vec![
            [0.,0.,5.],[1.,0.,6.],[0.,1.,4.],[-1.,0.,6.],[0.,-1.,4.],
        ],triangles:vec![[0,1,2],[0,2,3],[0,3,4],[0,4,1]]};
        let output=generate_contours(&s,0.5,5,0.,1e-8).unwrap();
        assert!(!output.is_empty());
        assert!(output.iter().all(|c|c.chain.vertices.len()>=2));
    }
    #[test]
    fn rejects_bad_intervals() {
        assert!(generate_contours(&ramp(),0.,5,0.,1e-6).is_err());
        assert!(generate_contours(&ramp(),0.5,0,0.,1e-6).is_err());
        assert!(generate_contours(&ramp(),0.5,5,0.,0.).is_err());
    }
    #[test]
    fn empty_surface_is_ok() {
        let out=generate_contours(&Surface{nodes:vec![],triangles:vec![]},0.5,5,0.,1e-6).unwrap();
        assert!(out.is_empty());
    }
    #[test]
    fn large_level_count_is_rejected() {
        assert!(generate_contours(&ramp(),1e-6,5,0.,1e-6).is_err());
    }
}
