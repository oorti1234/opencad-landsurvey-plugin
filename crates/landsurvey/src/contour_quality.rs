//! Read-only quality checks for assembled contour geometry.
use crate::contour_chains::ContourChain;
use std::collections::BTreeMap;

#[derive(Debug,Clone,PartialEq)]
pub struct ContourQuality {
    pub chains: usize,
    pub open_chains: usize,
    pub closed_chains: usize,
    pub edges: usize,
    pub zero_length_edges: usize,
    pub duplicate_edges: usize,
    pub invalid_chains: usize,
    pub malformed_closed_chains: usize,
}

/// Inspect contour chains without modifying or reordering their vertices.
/// Coincident edges are checked within each elevation separately.
pub fn audit_chains(chains:&[ContourChain], tolerance:f64)->Result<ContourQuality,String>{
    if !tolerance.is_finite() || tolerance<=0.0 {
        return Err("invalid contour audit XY tolerance".into());
    }
    let mut result=ContourQuality {
        chains:chains.len(),open_chains:0,closed_chains:0,edges:0,
        zero_length_edges:0,duplicate_edges:0,invalid_chains:0,malformed_closed_chains:0
    };
    let mut seen:BTreeMap<(u64,(i64,i64),(i64,i64)),usize>=BTreeMap::new();
    let key=|xy:[f64;2]|->Option<(i64,i64)> {
        if xy.iter().any(|n|!n.is_finite()){return None;}
        let x=(xy[0]/tolerance).round();let y=(xy[1]/tolerance).round();
        if !x.is_finite()||!y.is_finite()||x<i64::MIN as f64||x>=i64::MAX as f64
            ||y<i64::MIN as f64||y>=i64::MAX as f64{return None;}
        Some((x as i64,y as i64))
    };
    for chain in chains {
        if chain.closed {result.closed_chains+=1;}else{result.open_chains+=1;}
        if !chain.elevation.is_finite()||chain.vertices.len()<2 {
            result.invalid_chains+=1;continue;
        }
        let vertices:Option<Vec<_>>=chain.vertices.iter().copied().map(key).collect();
        let Some(vertices)=vertices else {result.invalid_chains+=1;continue;};
        let is_duplicate_end=vertices.first()==vertices.last();
        if chain.closed && (!is_duplicate_end || vertices.len()<4) {
            result.malformed_closed_chains+=1;
        }
        let level=chain.elevation.to_bits();
        for pair in vertices.windows(2) {
            result.edges+=1;
            if pair[0]==pair[1] {result.zero_length_edges+=1;continue;}
            let (a,b)=if pair[0]<pair[1]{(pair[0],pair[1])}else{(pair[1],pair[0])};
            if seen.insert((level,a,b),1).is_some(){result.duplicate_edges+=1;}
        }
        // A CAD closed polyline may store its final closure implicitly.
        if chain.closed && !is_duplicate_end {
            result.edges+=1;
            let a=*vertices.last().unwrap();let b=vertices[0];
            if a==b{result.zero_length_edges+=1;}
            else {
                let edge=if a<b{(level,a,b)}else{(level,b,a)};
                if seen.insert(edge,1).is_some(){result.duplicate_edges+=1;}
            }
        }
    }
    Ok(result)
}

/// Classify open-chain endpoints against actual TIN boundary edges and
/// previously identified ambiguous junctions. No geometry is changed.
#[derive(Debug,Clone,PartialEq)]
pub struct EndpointAudit {
    pub boundary:usize,
    pub junction:usize,
    pub interior:usize,
    pub invalid:usize,
    pub interior_samples:Vec<[f64;3]>,
}
pub fn audit_endpoints(
    surface:&crate::surface::Surface,
    chains:&[ContourChain],
    junctions:&[crate::contour_topology::Junction],
    tolerance:f64,
)->Result<EndpointAudit,String>{
    if !tolerance.is_finite() || tolerance<=0.0 {return Err("invalid endpoint tolerance".into());}
    let mut edges:BTreeMap<(usize,usize),usize>=BTreeMap::new();
    for tri in &surface.triangles {
        if tri.iter().any(|&i|i>=surface.nodes.len()) {return Err("invalid TIN triangle".into());}
        for (a,b) in [(tri[0],tri[1]),(tri[1],tri[2]),(tri[2],tri[0])] {
            let edge=if a<b{(a,b)}else{(b,a)};
            *edges.entry(edge).or_default()+=1;
        }
    }
    let boundary:Vec<_>=edges.into_iter().filter(|(_,n)|*n==1)
        .map(|((a,b),_)|(surface.nodes[a],surface.nodes[b])).collect();
    let mut out=EndpointAudit{boundary:0,junction:0,interior:0,invalid:0,interior_samples:vec![]};
    for c in chains.iter().filter(|c|!c.closed) {
        for &p in [c.vertices.first(),c.vertices.last()].into_iter().flatten() {
            if !c.elevation.is_finite() || p.iter().any(|v|!v.is_finite()) {
                out.invalid+=1;continue;
            }
            let on_edge=boundary.iter().any(|(a,b)|{
                let vx=b[0]-a[0];let vy=b[1]-a[1];
                let len2=vx*vx+vy*vy;
                if !len2.is_finite() || len2<=0.0 {return false;}
                let t=(((p[0]-a[0])*vx+(p[1]-a[1])*vy)/len2).clamp(0.,1.);
                let x=a[0]+t*vx;let y=a[1]+t*vy;
                (p[0]-x).hypot(p[1]-y)<=tolerance
            });
            if on_edge {out.boundary+=1;continue;}
            if junctions.iter().any(|j|j.elevation==c.elevation
                && (p[0]-j.xy[0]).hypot(p[1]-j.xy[1])<=tolerance) {
                out.junction+=1;continue;
            }
            out.interior+=1;
            if out.interior_samples.len()<8 {
                out.interior_samples.push([p[0],p[1],c.elevation]);
            }
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn chain(points:Vec<[f64;2]>,closed:bool)->ContourChain{
        ContourChain{elevation:41.5,vertices:points,closed}
    }
    #[test]
    fn identifies_boundary_and_interior_endpoints() {
        let surface=crate::surface::Surface{
            nodes:vec![[0.,0.,0.],[4.,0.,0.],[0.,4.,0.]],
            triangles:vec![[0,1,2]],
        };
        let chains=vec![chain(vec![[2.,0.],[1.,1.]],false)];
        let report=audit_endpoints(&surface,&chains,&[],1e-6).unwrap();
        assert_eq!(report.boundary,1);
        assert_eq!(report.interior,1);
    }
    #[test]
    fn identifies_ambiguous_interior_junction() {
        let surface=crate::surface::Surface{
            nodes:vec![[0.,0.,0.],[4.,0.,0.],[0.,4.,0.]],
            triangles:vec![[0,1,2]],
        };
        let j=crate::contour_topology::Junction{elevation:41.5,xy:[1.,1.],degree:4};
        let report=audit_endpoints(&surface,&[chain(vec![[2.,0.],[1.,1.]],false)],&[j],1e-6).unwrap();
        assert_eq!(report.boundary,1);
        assert_eq!(report.junction,1);
        assert_eq!(report.interior,0);
    }
    #[test]
    fn accepts_one_open_and_one_closed_contour(){
        let curves=vec![
            chain(vec![[0.,0.],[1.,0.]],false),
            chain(vec![[2.,0.],[3.,0.],[3.,1.],[2.,0.]],true)
        ];
        let report=audit_chains(&curves,1e-6).unwrap();
        assert_eq!(report.edges,4);
        assert_eq!(report.open_chains,1);
        assert_eq!(report.closed_chains,1);
        assert_eq!(report.duplicate_edges,0);
        assert_eq!(report.malformed_closed_chains,0);
    }
    #[test]
    fn detects_duplicate_and_collapsed_segments(){
        let curves=vec![
            chain(vec![[0.,0.],[1.,0.],[1.,0.]],false),
            chain(vec![[1.,0.],[0.,0.]],false)
        ];
        let report=audit_chains(&curves,1e-6).unwrap();
        assert_eq!(report.duplicate_edges,1);
        assert_eq!(report.zero_length_edges,1);
    }
    #[test]
    fn detects_invalid_closed_topology_without_mutation(){
        let curves=vec![chain(vec![[0.,0.],[1.,0.],[1.,1.]],true)];
        let original=curves.clone();
        let report=audit_chains(&curves,1e-6).unwrap();
        assert_eq!(report.malformed_closed_chains,1);
        assert_eq!(curves,original);
    }
}
