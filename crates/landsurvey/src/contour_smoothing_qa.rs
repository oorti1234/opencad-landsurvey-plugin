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
/// Directed maximum distance of polyline vertices to the opposite polyline.
/// For two straight polylines this gives the exact directed Hausdorff distance:
/// distance to a union of segments is 1-Lipschitz, but its maximum along an
/// individual segment can occur between endpoints. Sample each source
/// segment at evenly spaced fractions to expose deviations at cut corners.
/// This is a sampled estimate, not a certified Hausdorff bound.
fn sampled_directed_distance(a:&[[f64;2]],b:&[[f64;2]])->f64 {
    let mut worst:f64=0.;
    for edge in a.windows(2) {
        for k in 0..=8 {
            let t=k as f64/8.;
            let p=[edge[0][0]*(1.-t)+edge[1][0]*t,
                   edge[0][1]*(1.-t)+edge[1][1]*t];
            let d=b.windows(2).map(|other|
                point_to_segment(p,other[0],other[1])
            ).fold(f64::INFINITY,f64::min);
            worst=worst.max(d);
        }
    }
    worst
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
        let original_to_preview=sampled_directed_distance(
            &original.vertices,&smoothed.vertices
        );
        let preview_to_original=sampled_directed_distance(
            &smoothed.vertices,&original.vertices
        );
        max_distance=max_distance.max(original_to_preview).max(preview_to_original);
        vertices+=smoothed.vertices.len();
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
/// Intersection audit between contours at different elevations.
/// Detects crossings, endpoint touches, and positive-length collinear overlaps.
/// A tolerance is expressed in XY drawing units; near-misses are not classified.
#[derive(Debug,Clone,PartialEq,Default)]
pub struct InterLevelConflicts {
    pub crossings:usize,
    pub touches:usize,
    pub overlaps:usize,
}
fn orient(a:[f64;2],b:[f64;2],p:[f64;2])->f64 {
    (b[0]-a[0])*(p[1]-a[1])-(b[1]-a[1])*(p[0]-a[0])
}
fn on_segment(p:[f64;2],a:[f64;2],b:[f64;2],tol:f64)->bool {
    point_to_segment(p,a,b)<=tol
}
pub fn inter_level_conflicts(
    curves:&[ContourChain],tol:f64
)->Result<InterLevelConflicts,String> {
    if !tol.is_finite()||tol<=0. {return Err("invalid XY tolerance".into());}
    let mut out=InterLevelConflicts::default();
    for c in curves {
        if !c.elevation.is_finite()||c.vertices.len()<2
            ||c.vertices.iter().flatten().any(|v|!v.is_finite()){
            return Err("invalid contour geometry".into());
        }
    }
    for i in 0..curves.len(){
        for j in i+1..curves.len(){
            if curves[i].elevation==curves[j].elevation {continue;}
            for ab in curves[i].vertices.windows(2){
                for cd in curves[j].vertices.windows(2){
                    let (a,b,c,d)=(ab[0],ab[1],cd[0],cd[1]);
                    let len_ab=(b[0]-a[0]).hypot(b[1]-a[1]);
                    let len_cd=(d[0]-c[0]).hypot(d[1]-c[1]);
                    if len_ab<=tol||len_cd<=tol {continue;}
                    let cross1=orient(a,b,c);let cross2=orient(a,b,d);
                    let cross3=orient(c,d,a);let cross4=orient(c,d,b);
                    let collinear=cross1.abs()<=tol*len_ab
                        &&cross2.abs()<=tol*len_ab
                        &&cross3.abs()<=tol*len_cd
                        &&cross4.abs()<=tol*len_cd;
                    if collinear {
                        let (u,v,w,x)=if (b[0]-a[0]).abs()>=(b[1]-a[1]).abs(){
                            (a[0],b[0],c[0],d[0])
                        } else {(a[1],b[1],c[1],d[1])};
                        let overlap=u.max(v).min(w.max(x))-u.min(v).max(w.min(x));
                        if overlap>tol {out.overlaps+=1;}
                        else if overlap>=-tol {out.touches+=1;}
                    } else if proper_cross(a,b,c,d) {
                        out.crossings+=1;
                    } else if [a,b].iter().any(|&p|on_segment(p,c,d,tol))
                        ||[c,d].iter().any(|&p|on_segment(p,a,b,tol)) {
                        out.touches+=1;
                    }
                }
            }
        }
    }
    Ok(out)
}
/// Sample each preview polyline against the TIN at vertices and segment
/// midpoints. Fail closed when a sample lies outside the triangulation.
#[derive(Debug,Clone,PartialEq)]
pub struct TinElevationQuality {
    pub samples:usize,
    pub outside_tin:usize,
    pub max_absolute_z_error:f64,
}
pub fn audit_tin_elevation(
    surface:&crate::surface::Surface,
    previews:&[ContourChain],
)->Result<TinElevationQuality,String>{
    if surface.triangles.iter().flatten().any(|&i|i>=surface.nodes.len()) {
        return Err("invalid TIN indices".into());
    }
    let index=crate::tin_spatial_index::TinGridIndex::build(surface)?;
    let mut result=TinElevationQuality{samples:0,outside_tin:0,max_absolute_z_error:0.};
    for chain in previews {
        if !chain.elevation.is_finite() || chain.vertices.len()<2 ||
           chain.vertices.iter().flatten().any(|v|!v.is_finite()) {
            return Err("invalid contour samples".into());
        }
        let mut samples=chain.vertices.clone();
        samples.extend(chain.vertices.windows(2).map(|v|
            [(v[0][0]+v[1][0])*0.5,(v[0][1]+v[1][1])*0.5]));
        for p in samples {
            result.samples+=1;
            match index.interpolate_z(surface,p[0],p[1]).or_else(||surface.interpolate_z(p[0],p[1])) {
                Some(z) if z.is_finite()=>{
                    result.max_absolute_z_error=result.max_absolute_z_error.max(
                        (z-chain.elevation).abs()
                    );
                }
                _=>result.outside_tin+=1,
            }
        }
    }
    Ok(result)
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
    fn chaikin_corner_cut_has_nonzero_bidirectional_distance(){
        let source=[c(1.,vec![[0.,0.],[1.,1.],[2.,0.]])];
        let preview=[crate::contour_smoothing::smooth_preview(&source[0],1).unwrap()];
        let q=audit_smoothing(&source,&preview).unwrap();
        assert!(q.max_vertex_displacement > 0.1, "deviation: {}", q.max_vertex_displacement);
    }
    #[test]
    fn detects_deviation_between_vertices(){
        let source=[c(1.,vec![[0.,0.],[1.,1.],[2.,0.]])];
        let preview=[c(1.,vec![[0.,0.],[2.,0.]])];
        let q=audit_smoothing(&source,&preview).unwrap();
        assert!(q.max_vertex_displacement>0.9);
    }
    #[test]
    fn flags_crossing_between_elevations(){
        let src=[c(1.,vec![[0.,0.],[2.,2.]]),c(2.,vec![[0.,2.],[2.,0.]])];
        let q=audit_smoothing(&src,&src).unwrap();
        assert_eq!(q.crossings_between_levels,1);
    }
    #[test]
    fn finds_off_contour_vertical_deviation(){
        let tin=crate::surface::Surface {
            nodes:vec![[0.,0.,0.],[2.,0.,2.],[0.,2.,0.]],
            triangles:vec![[0,1,2]]
        };
        let curve=[c(1.,vec![[0.,0.],[1.,0.]])];
        let report=audit_tin_elevation(&tin,&curve).unwrap();
        assert!((report.max_absolute_z_error-1.).abs()<1e-9);
        assert_eq!(report.outside_tin,0);
    }
    #[test]
    fn reports_outside_tin_samples(){
        let tin=crate::surface::Surface {
            nodes:vec![[0.,0.,0.],[2.,0.,2.],[0.,2.,0.]],
            triangles:vec![[0,1,2]]
        };
        let curve=[c(1.,vec![[3.,0.],[4.,0.]])];
        assert_eq!(audit_tin_elevation(&tin,&curve).unwrap().outside_tin,3);
    }
    #[test]
    fn finds_touch_overlap_and_crossing_between_levels(){
        let touch=[c(1.,vec![[0.,0.],[2.,0.]]),
                   c(2.,vec![[2.,0.],[2.,2.]])];
        let overlap=[c(1.,vec![[0.,0.],[3.,0.]]),
                     c(2.,vec![[1.,0.],[2.,0.]])];
        let crossing=[c(1.,vec![[0.,0.],[2.,2.]]),
                      c(2.,vec![[0.,2.],[2.,0.]])];
        assert_eq!(inter_level_conflicts(&touch,1e-8).unwrap().touches,1);
        assert_eq!(inter_level_conflicts(&overlap,1e-8).unwrap().overlaps,1);
        assert_eq!(inter_level_conflicts(&crossing,1e-8).unwrap().crossings,1);
    }
    #[test]
    fn refuses_elevation_changes(){
        let a=[c(1.,vec![[0.,0.],[2.,0.]])];
        let b=[c(2.,vec![[0.,0.],[2.,0.]])];
        assert!(audit_smoothing(&a,&b).is_err());
    }
}
