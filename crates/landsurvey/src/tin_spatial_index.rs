//! Reusable uniform XY grid for TIN interpolation. Candidate triangles are
//! visited in original index order to preserve the first-hit behavior of
//! Surface::interpolate_z, including points on shared triangle boundaries.
use crate::surface::Surface;
#[derive(Debug,Clone)]
pub struct TinGridIndex {
    min:[f64;2],cell:[f64;2],nx:usize,ny:usize,
    cells:Vec<Vec<usize>>,
}
impl TinGridIndex {
    pub fn build(surface:&Surface)->Result<Self,String>{
        if surface.nodes.is_empty() || surface.triangles.is_empty(){
            return Err("TIN is empty".into());
        }
        if surface.nodes.iter().flatten().any(|x|!x.is_finite())
            ||surface.triangles.iter().flatten().any(|&i|i>=surface.nodes.len()){
            return Err("invalid TIN coordinates or indices".into());
        }
        let (mut minx,mut maxx,mut miny,mut maxy)=(
            f64::INFINITY,f64::NEG_INFINITY,f64::INFINITY,f64::NEG_INFINITY
        );
        for p in &surface.nodes{
            minx=minx.min(p[0]);maxx=maxx.max(p[0]);
            miny=miny.min(p[1]);maxy=maxy.max(p[1]);
        }
        let n=((surface.triangles.len() as f64).sqrt().ceil() as usize).clamp(1,64);
        let nx=n;let ny=n;
        let cell=[((maxx-minx)/nx as f64).max(1e-9),
                  ((maxy-miny)/ny as f64).max(1e-9)];
        let mut cells=vec![Vec::new();nx*ny];
        let slot=|x:f64,y:f64|->(usize,usize) {
            let ix=(((x-minx)/cell[0]).floor() as isize).clamp(0,nx as isize-1) as usize;
            let iy=(((y-miny)/cell[1]).floor() as isize).clamp(0,ny as isize-1) as usize;
            (ix,iy)
        };
        for (idx,t) in surface.triangles.iter().enumerate(){
            let pts=[surface.nodes[t[0]],surface.nodes[t[1]],surface.nodes[t[2]]];
            let lx=pts.iter().map(|p|p[0]).fold(f64::INFINITY,f64::min);
            let hx=pts.iter().map(|p|p[0]).fold(f64::NEG_INFINITY,f64::max);
            let ly=pts.iter().map(|p|p[1]).fold(f64::INFINITY,f64::min);
            let hy=pts.iter().map(|p|p[1]).fold(f64::NEG_INFINITY,f64::max);
            let (ax,ay)=slot(lx,ly);let (bx,by)=slot(hx,hy);
            for iy in ay..=by {for ix in ax..=bx {
                cells[iy*nx+ix].push(idx);
            }}
        }
        Ok(Self{min:[minx,miny],cell,nx,ny,cells})
    }
    pub fn interpolate_z(&self,surface:&Surface,x:f64,y:f64)->Option<f64>{
        if !x.is_finite()||!y.is_finite(){return None;}
        let fx=(x-self.min[0])/self.cell[0];
        let fy=(y-self.min[1])/self.cell[1];
        if !fx.is_finite()||!fy.is_finite(){return None;}
        let ix=fx.floor() as isize;let iy=fy.floor() as isize;
        // Probe adjacent cells to preserve the original barycentric boundary
        // tolerance near grid edges.
        let mut candidates=std::collections::BTreeSet::new();
        for yy in iy.saturating_sub(1)..=iy.saturating_add(1) {
            for xx in ix.saturating_sub(1)..=ix.saturating_add(1) {
                if xx>=0&&yy>=0&&(xx as usize)<self.nx&&(yy as usize)<self.ny {
                    candidates.extend(self.cells[yy as usize*self.nx+xx as usize].iter().copied());
                }
            }
        }
        for idx in candidates{
            let t=surface.triangles[idx];
            let (a,b,c)=(surface.nodes[t[0]],surface.nodes[t[1]],surface.nodes[t[2]]);
            let d=(b[1]-c[1])*(a[0]-c[0])+(c[0]-b[0])*(a[1]-c[1]);
            if d.abs()<1e-12{continue;}
            let wa=((b[1]-c[1])*(x-c[0])+(c[0]-b[0])*(y-c[1]))/d;
            let wb=((c[1]-a[1])*(x-c[0])+(a[0]-c[0])*(y-c[1]))/d;
            let wc=1.-wa-wb;
            if wa>=-1e-9&&wb>=-1e-9&&wc>=-1e-9{
                return Some(wa*a[2]+wb*b[2]+wc*c[2]);
            }
        }
        None
    }
}
#[cfg(test)]
mod tests{
    use super::*;
    #[test]
    fn indexed_interpolation_matches_scan_on_grid_and_edges(){
        let surface=Surface{
            nodes:vec![[0.,0.,0.],[2.,0.,2.],[2.,2.,4.],[0.,2.,2.]],
            triangles:vec![[0,1,2],[0,2,3]]
        };
        let index=TinGridIndex::build(&surface).unwrap();
        for yi in 0..=20{for xi in 0..=20{
            let x=xi as f64/10.;let y=yi as f64/10.;
            let a=surface.interpolate_z(x,y);
            let b=index.interpolate_z(&surface,x,y);
            assert_eq!(a.is_some(),b.is_some());
            if let (Some(a),Some(b))=(a,b){assert!((a-b).abs()<1e-10);}
        }}
    }
    #[test]
    fn outside_and_invalid_samples_are_rejected(){
        let surface=Surface{nodes:vec![[0.,0.,0.],[1.,0.,1.],[0.,1.,1.]],
            triangles:vec![[0,1,2]]};
        let idx=TinGridIndex::build(&surface).unwrap();
        assert_eq!(idx.interpolate_z(&surface,3.,3.),None);
        assert_eq!(idx.interpolate_z(&surface,f64::NAN,0.),None);
    }
}
