//! Candidate pairings for four contour branches at an exact-level TIN saddle.
//! Both non-crossing interpretations are retained: there is no canonical
//! physical connection on the unperturbed piecewise-linear surface.
use crate::contours::ContourSegment;

#[derive(Debug, Clone, PartialEq)]
pub struct SaddleOptions {
    /// Two candidate pairings of the original input segment indices.
    pub alternatives: [[(usize, usize); 2]; 2],
}

/// Enumerate two locally non-crossing pairing choices for exactly four
/// independent branches incident at the same point/elevation. Does NOT
/// choose a pairing or alter CAD geometry.
pub fn saddle_pairing_options(
    segments: &[ContourSegment],
    center: [f64; 2],
    elevation: f64,
    tolerance: f64,
) -> Option<SaddleOptions> {
    if segments.len()!=4 || !elevation.is_finite()
        || !tolerance.is_finite() || tolerance<=0.0
        || center.iter().any(|v|!v.is_finite()) { return None; }
    let mut angles=Vec::with_capacity(4);
    for (index,s) in segments.iter().enumerate() {
        if s.elevation!=elevation {return None;}
        let a=(s.start[0]-center[0]).hypot(s.start[1]-center[1])<=tolerance;
        let b=(s.end[0]-center[0]).hypot(s.end[1]-center[1])<=tolerance;
        if a==b {return None;}
        let other=if a{s.end}else{s.start};
        let dx=other[0]-center[0]; let dy=other[1]-center[1];
        if !dx.is_finite() || !dy.is_finite() || dx.hypot(dy)<=tolerance {
            return None;
        }
        angles.push((dy.atan2(dx),index));
    }
    angles.sort_by(|a,b|a.0.total_cmp(&b.0));
    let i:Vec<usize>=angles.iter().map(|(_,index)|*index).collect();
    Some(SaddleOptions{alternatives:[
        [(i[0],i[1]),(i[2],i[3])],
        [(i[1],i[2]),(i[3],i[0])],
    ]})
}
/// Only offer the two candidate connections when the actual TIN confirms a
/// complete alternating interior saddle at this elevation. A confirmed saddle
/// still has two valid local resolutions; choosing one requires an explicit
/// policy, not an arbitrary angle-based guess.
pub fn verified_saddle_options(
    surface: &crate::surface::Surface,
    segments: &[ContourSegment],
    center: [f64; 2],
    elevation: f64,
    tolerance: f64,
) -> Option<SaddleOptions> {
    if !crate::contour_topology::is_exact_level_saddle(
        surface, center, elevation, tolerance
    ) {
        return None;
    }
    saddle_pairing_options(segments, center, elevation, tolerance)
}

/// Read-only geometric preview of a pairing: each proposed continuation
/// consists of two original branch endpoints joined through the same saddle.
/// The preview does not change the terrain or decide which option is correct.
#[derive(Debug, Clone, PartialEq)]
pub struct SaddlePreview {
    pub alternatives: [[[f64; 2]; 3]; 4],
}

/// Returns four 3-point preview paths: the first two belong to choice A
/// and the last two to choice B.
pub fn saddle_preview_paths(
    segments: &[ContourSegment],
    center: [f64; 2],
    elevation: f64,
    tolerance: f64,
) -> Option<SaddlePreview> {
    let choices = saddle_pairing_options(segments, center, elevation, tolerance)?;
    let endpoint = |index: usize| -> [f64; 2] {
        let s = &segments[index];
        if (s.start[0]-center[0]).hypot(s.start[1]-center[1]) <= tolerance {
            s.end
        } else { s.start }
    };
    let make = |a: usize, b: usize| [endpoint(a), center, endpoint(b)];
    let a = choices.alternatives[0];
    let b = choices.alternatives[1];
    Some(SaddlePreview { alternatives: [
        make(a[0].0,a[0].1), make(a[1].0,a[1].1),
        make(b[0].0,b[0].1), make(b[1].0,b[1].1),
    ]})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn four_branches_have_two_distinct_non_crossing_options() {
        let c=[355420.383,928758.559];
        let s=vec![
            ContourSegment{elevation:41.5,start:c,end:[c[0]+1.,c[1]]},
            ContourSegment{elevation:41.5,start:c,end:[c[0],c[1]+1.]},
            ContourSegment{elevation:41.5,start:c,end:[c[0]-1.,c[1]]},
            ContourSegment{elevation:41.5,start:c,end:[c[0],c[1]-1.]},
        ];
        let o=saddle_pairing_options(&s,c,41.5,1e-6).unwrap();
        assert_ne!(o.alternatives[0],o.alternatives[1]);
        assert_eq!(o.alternatives[0],[(3,0),(1,2)]);
        assert_eq!(o.alternatives[1],[(0,1),(2,3)]);
    }
    #[test]
    fn verified_options_require_real_tin_saddle() {
        let nodes=vec![
            [0.,0.,5.],[1.,0.,6.],[0.,1.,4.],[-1.,0.,6.],[0.,-1.,4.]
        ];
        let closed=crate::surface::Surface {
            nodes:nodes.clone(),
            triangles:vec![[0,1,2],[0,2,3],[0,3,4],[0,4,1]]
        };
        let open=crate::surface::Surface {
            nodes,
            triangles:vec![[0,1,2],[0,2,3],[0,3,4]]
        };
        let segments=vec![
            ContourSegment{elevation:5.,start:[0.,0.],end:[0.5,0.5]},
            ContourSegment{elevation:5.,start:[0.,0.],end:[-0.5,0.5]},
            ContourSegment{elevation:5.,start:[0.,0.],end:[-0.5,-0.5]},
            ContourSegment{elevation:5.,start:[0.,0.],end:[0.5,-0.5]},
        ];
        assert!(verified_saddle_options(&closed,&segments,[0.,0.],5.,1e-6).is_some());
        assert!(verified_saddle_options(&open,&segments,[0.,0.],5.,1e-6).is_none());
        assert!(verified_saddle_options(&closed,&segments,[0.,0.],5.5,1e-6).is_none());
    }
    #[test]
    fn preview_paths_preserve_original_endpoints_and_center() {
        let center=[0.,0.];
        let branches=[
            ContourSegment{elevation:5.,start:center,end:[1.,0.]},
            ContourSegment{elevation:5.,start:center,end:[0.,1.]},
            ContourSegment{elevation:5.,start:center,end:[-1.,0.]},
            ContourSegment{elevation:5.,start:center,end:[0.,-1.]},
        ];
        let preview=saddle_preview_paths(&branches,center,5.,1e-6).unwrap();
        assert!(preview.alternatives.iter().all(|p|p[1]==center));
        let mut endpoints=preview.alternatives[0..2].iter()
            .flat_map(|p|[p[0],p[2]]).collect::<Vec<_>>();
        endpoints.sort_by(|a,b| a[0].total_cmp(&b[0]).then(a[1].total_cmp(&b[1])));
        assert_eq!(endpoints,vec![[-1.,0.],[0.,-1.],[0.,1.],[1.,0.]]);
    }
    #[test]
    fn rejects_non_saddle_or_wrong_level() {
        let s=vec![ContourSegment{elevation:41.5,start:[0.,0.],end:[1.,0.]};4];
        assert!(saddle_pairing_options(&s,[0.,0.],40.,1e-6).is_none());
        assert!(saddle_pairing_options(&s,[0.,0.],41.5,0.).is_none());
    }
}
