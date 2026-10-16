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
    fn rejects_non_saddle_or_wrong_level() {
        let s=vec![ContourSegment{elevation:41.5,start:[0.,0.],end:[1.,0.]};4];
        assert!(saddle_pairing_options(&s,[0.,0.],40.,1e-6).is_none());
        assert!(saddle_pairing_options(&s,[0.,0.],41.5,0.).is_none());
    }
}
