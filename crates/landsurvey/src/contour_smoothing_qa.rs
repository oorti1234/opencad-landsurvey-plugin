//! Read-only comparison of cartographic smoothing against source contour lines.
use crate::contour_chains::ContourChain;
#[derive(Debug,Clone,PartialEq)]
pub struct SmoothingQuality {
    pub max_vertex_displacement:f64,
    pub crossings_between_levels:usize,
    pub compared_vertices:usize,
}
fn point_to_segment(p:[f64;2],a:[f64;2],b:[f64;2])->f64{
    let dx=b[0]-a[0];let dy=b[1]-a[1];let len=dx*dx+dy*dy;
    if len==0. {return (p[0]-a[0]).hypot(p[1]-a[1]);}
    let t=(((p[0]-a[0])*dx+(p[1]-a[1])*dy)/len).clamp(0.,1.);
    (p[0]-a[0]-t*dx).hypot(p[1]-a[1]-t*dy)
}
fn cross(a:[f64;2],b:[f64;2],c:[f64;2])->f64{
    (b[0]-a[0])*(c[1]-a[1])-(b[1]-a[1])*(c[0]-a[0])
}
fn proper_cross(a:[f64;2],b:[f64;2],c:[f64;2],d:[f64;2])->bool{
    let x=cross(a,b,c);let y=cross(a,b,d);
    let z=cross(c,d,a);let w=cross(c,d,b);
    (x>0. && y<0. || x<0. && y>0.) &&
    (z>0. && w<0. || z<0. && w>0.)
}
/// Compare a preview with the corresponding source chains. Crossing detection
/// counts strict segment intersections across differing elevations; tangencies
/// and collinear overlaps require additional checks before certification.
pub fn audit_smoothing(source:&[ContourChain],preview:&[ContourChain])->Result<SmoothingQuality,String>{
    if source.len()!=preview.len(){return Err("source and preview chain counts differ".into());}
    let mut max_distance:f64=0.;
    let mut vertices=0usize;
    for (original,smoothed) in source.iter().zip(preview){
        if original.elevation!=smoothed.elevation || original.closed!=smoothed.closed {
            return Err("preview changes elevation or closure".into());
        }
        if original.vertices.len()<2 || smoothed.vertices.len()<2 ||
            original.vertices.iter().chain(smoothed.vertices.iter())
                .flatten().any(|n|!n.is_finite()) {
            return Err("invalid contour geometry".into());
        }
        for &p in &smoothed.vertices {
            let dist=original.vertices.windows(2)
                .map(|segment|point_to_segment(p,segment[0],segment[1]))
                .fold(f64::INFINITY,f64::min);
            max_distance=max_distance.max(dist);
            vertices+=1;
        }
    }
    let mut crossings=0usize;
    for i in 0..preview.len(){
        for j in i+1..preview.len(){
            if preview[i].elevation==preview[j].elevation {continue;}
            for a in preview[i].vertices.windows(2){
                for b in preview[j].vertices.windows(2){
                    if proper_cross(a[0],a[1],b[0],b[1]){crossings+=1;}
                }
            }
        }
    }
    Ok(SmoothingQuality{max_vertex_displacement:max_distance,
        crossings_between_levels:crossings,compared_vertices:vertices})
}
#[cfg(test)]
mod tests{
    use super::*;
    fn c(z:f64,p:Vec<[f64;2]>)->ContourChain{ContourChain{elevation:z,vertices:p,closed:false}}
    #[test]
    fn reports_displacement(){
        let source=[c(1.,vec![[0.,0.],[2.,0.]])];
        let preview=[c(1.,vec![[0.,0.],[1.,0.25],[2.,0.]])];
        let q=audit_smoothing(&source,&preview).unwrap();
        assert!((q.max_vertex_displacement-0.25).abs()<1e-10);
        assert_eq!(q.crossings_between_levels,0);
    }
    #[test]
    fn flags_crossing_between_elevations(){
        let src=[c(1.,vec![[0.,0.],[2.,2.]]),c(2.,vec![[0.,2.],[2.,0.]])];
        let q=audit_smoothing(&src,&src).unwrap();
        assert_eq!(q.crossings_between_levels,1);
    }
    #[test]
    fn refuses_elevation_changes(){
        let a=[c(1.,vec![[0.,0.],[2.,0.]])];
        let b=[c(2.,vec![[0.,0.],[2.,0.]])];
        assert!(audit_smoothing(&a,&b).is_err());
    }
}
