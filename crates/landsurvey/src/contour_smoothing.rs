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
    fn rejects_excessive_passes(){
        let c=ContourChain{elevation:42.,vertices:vec![[0.,0.],[1.,1.]],closed:false};
        assert!(smooth_preview(&c,4).is_err());
    }
}
