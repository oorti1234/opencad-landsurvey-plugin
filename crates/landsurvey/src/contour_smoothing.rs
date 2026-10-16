//! Non-destructive, corner-cutting previews of TIN contour chains.
//! Smoothed coordinates are visualization only, never survey-grade isolines.
use crate::contour_chains::ContourChain;

/// One Chaikin corner-cutting pass. Open chain endpoints remain fixed;
/// closed chains retain closure. Never changes original input.
pub fn smooth_preview(chain:&ContourChain,passes:usize)->Result<ContourChain,String>{
    if passes>3 {return Err("maximum three smoothing passes".into());}
    if !chain.elevation.is_finite() || chain.vertices.iter().flatten().any(|v|!v.is_finite()){
        return Err("non-finite contour geometry".into());
    }
    if chain.vertices.len()<2 {return Err("contour has fewer than two vertices".into());}
    let mut points=chain.vertices.clone();
    for _ in 0..passes {
        if chain.closed {
            let ring=&points[..if points.first()==points.last(){points.len()-1}else{points.len()}];
            if ring.len()<3{return Err("closed contour has too few vertices".into());}
            let mut next=Vec::with_capacity(ring.len()*2+1);
            for i in 0..ring.len(){
                let a=ring[i];let b=ring[(i+1)%ring.len()];
                next.push([0.75*a[0]+0.25*b[0],0.75*a[1]+0.25*b[1]]);
                next.push([0.25*a[0]+0.75*b[0],0.25*a[1]+0.75*b[1]]);
            }
            next.push(next[0]);
            points=next;
        } else {
            let mut next=Vec::with_capacity(points.len()*2);
            next.push(points[0]);
            for pair in points.windows(2){
                let a=pair[0];let b=pair[1];
                next.push([0.75*a[0]+0.25*b[0],0.75*a[1]+0.25*b[1]]);
                next.push([0.25*a[0]+0.75*b[0],0.25*a[1]+0.75*b[1]]);
            }
            next.push(*points.last().unwrap());
            points=next;
        }
        if points.len()>200_000{return Err("smoothing preview exceeds vertex limit".into());}
    }
    Ok(ContourChain{elevation:chain.elevation,vertices:points,closed:chain.closed})
}
/// First adaptive stage: accept smoothing independently for each complete
/// contour; retain its exact original coordinates if the sampled XY deviation
/// exceeds the requested limit. Never changes topology or elevation.
pub fn adaptive_preview(
    chains:&[ContourChain],
    passes:usize,
    max_xy:f64
)->Result<(Vec<ContourChain>,usize,usize),String>{
    if !max_xy.is_finite() || max_xy<=0. {return Err("invalid XY limit".into());}
    let mut output=Vec::with_capacity(chains.len());
    let mut smoothed_count=0;
    let mut unchanged_count=0;
    for chain in chains {
        let candidate=smooth_preview(chain,passes)?;
        let report=crate::contour_smoothing_qa::audit_smoothing(
            std::slice::from_ref(chain),std::slice::from_ref(&candidate)
        )?;
        if report.max_vertex_displacement<=max_xy {
            output.push(candidate);smoothed_count+=1;
        } else {
            output.push(chain.clone());unchanged_count+=1;
        }
    }
    Ok((output,smoothed_count,unchanged_count))
}

/// Experimental local corner cutting on open chains only. Each vertex has a
/// separately limited corner-cut fraction based on adjacent segment lengths.
/// Closed loops are retained exactly in this first phase to avoid seams.
pub fn local_preview(
    chains:&[ContourChain],max_xy:f64
)->Result<(Vec<ContourChain>,usize,usize),String>{
    if !max_xy.is_finite() || max_xy<=0. {return Err("invalid local XY limit".into());}
    let mut output=Vec::with_capacity(chains.len());
    let mut changed=0;
    let mut kept=0;
    for chain in chains {
        if chain.closed || chain.vertices.len()<3 {
            output.push(chain.clone());kept+=1;continue;
        }
        if !chain.elevation.is_finite() || chain.vertices.iter().flatten().any(|v|!v.is_finite()){
            return Err("invalid contour geometry".into());
        }
        let pts=&chain.vertices;
        let mut fractions=vec![0.;pts.len()];
        for i in 1..pts.len()-1 {
            let left=(pts[i][0]-pts[i-1][0]).hypot(pts[i][1]-pts[i-1][1]);
            let right=(pts[i+1][0]-pts[i][0]).hypot(pts[i+1][1]-pts[i][1]);
            let span=left.max(right);
            if span>0. {
                fractions[i]=0.25_f64.min(max_xy/(2.*span));
            }
        }
        let mut vertices=vec![pts[0]];
        for i in 0..pts.len()-1 {
            let a=pts[i];let b=pts[i+1];
            let fa=fractions[i];let fb=fractions[i+1];
            if fa>0. {vertices.push([a[0]*(1.-fa)+b[0]*fa,a[1]*(1.-fa)+b[1]*fa]);}
            if fb>0. {vertices.push([a[0]*fb+b[0]*(1.-fb),a[1]*fb+b[1]*(1.-fb)]);}
            if fb==0. {vertices.push(b);}
        }
        let candidate=ContourChain{elevation:chain.elevation,vertices,closed:false};
        let qa=crate::contour_smoothing_qa::audit_smoothing(
            std::slice::from_ref(chain),std::slice::from_ref(&candidate)
        )?;
        if qa.max_vertex_displacement<=max_xy {
            if candidate.vertices!=chain.vertices {changed+=1;} else {kept+=1;}
            output.push(candidate);
        } else {
            kept+=1;output.push(chain.clone());
        }
    }
    Ok((output,changed,kept))
}

#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn preserves_open_endpoints_and_height(){
        let c=ContourChain{elevation:41.5,vertices:vec![[0.,0.],[1.,1.],[2.,0.]],closed:false};
        let s=smooth_preview(&c,1).unwrap();
        assert_eq!(s.vertices.first(),c.vertices.first());
        assert_eq!(s.vertices.last(),c.vertices.last());
        assert_eq!(s.elevation,41.5);
        assert_eq!(c.vertices.len(),3);
    }
    #[test]
    fn closed_ring_stays_closed(){
        let c=ContourChain{elevation:42.,vertices:vec![[0.,0.],[2.,0.],[0.,2.],[0.,0.]],closed:true};
        let s=smooth_preview(&c,2).unwrap();
        assert!(s.closed);
        assert_eq!(s.vertices.first(),s.vertices.last());
    }
    #[test]
    fn adaptive_keeps_sharp_curve_when_limit_is_tiny() {
        let original=ContourChain{elevation:41.5,
            vertices:vec![[0.,0.],[1.,1.],[2.,0.]],closed:false};
        let (result,smooth,kept)=adaptive_preview(
            std::slice::from_ref(&original),1,0.001
        ).unwrap();
        assert_eq!((smooth,kept),(0,1));
        assert_eq!(result[0],original);
    }
    #[test]
    fn adaptive_smooths_when_limit_is_generous(){
        let original=ContourChain{elevation:41.5,
            vertices:vec![[0.,0.],[1.,1.],[2.,0.]],closed:false};
        let (result,smooth,kept)=adaptive_preview(
            std::slice::from_ref(&original),1,10.
        ).unwrap();
        assert_eq!((smooth,kept),(1,0));
        assert_ne!(result[0].vertices,original.vertices);
        assert_eq!(result[0].elevation,original.elevation);
    }
    #[test]
    fn local_smooths_sharp_corner_with_small_limit(){
        let c=ContourChain{elevation:41.5,
            vertices:vec![[0.,0.],[10.,10.],[20.,0.]],closed:false};
        let (out,changed,kept)=local_preview(&[c.clone()],0.5).unwrap();
        assert_eq!((changed,kept),(1,0));
        assert_eq!(out[0].vertices.first(),c.vertices.first());
        assert_eq!(out[0].vertices.last(),c.vertices.last());
        let qa=crate::contour_smoothing_qa::audit_smoothing(&[c],&out).unwrap();
        assert!(qa.max_vertex_displacement<=0.5);
    }
    #[test]
    fn local_leaves_closed_saddle_loops_untouched(){
        let c=ContourChain{elevation:5.,closed:true,
            vertices:vec![[0.,0.],[1.,0.],[1.,1.],[0.,0.]]};
        let (out,changed,kept)=local_preview(&[c.clone()],0.5).unwrap();
        assert_eq!((changed,kept),(0,1));
        assert_eq!(out[0],c);
    }
    #[test]
    fn rejects_excessive_passes(){
        let c=ContourChain{elevation:42.,vertices:vec![[0.,0.],[1.,1.]],closed:false};
        assert!(smooth_preview(&c,4).is_err());
    }
}
