//! Assemble independent isoline segments into ordered chains.
//! This module is experimental; CAD/host integration is deliberately separate.
use crate::contours::ContourSegment;
use std::collections::{BTreeMap, BTreeSet};

#[derive(Debug, Clone, PartialEq)]
pub struct ContourChain {
    pub elevation: f64,
    pub vertices: Vec<[f64; 2]>,
    pub closed: bool,
}

fn key(p: [f64; 2], tol: f64) -> Option<(i64, i64)> {
    if p.iter().any(|v| !v.is_finite()) { return None; }
    let x = (p[0] / tol).round();
    let y = (p[1] / tol).round();
    if x < i64::MIN as f64 || x > i64::MAX as f64 ||
       y < i64::MIN as f64 || y > i64::MAX as f64 { return None; }
    Some((x as i64, y as i64))
}

/// Join raw segments by elevation and endpoint coordinates. Returns an error
/// rather than guessing connectivity at branching or ambiguous nodes.
/// The caller chooses an XY tolerance expressed in drawing units.
pub fn stitch_segments(
    segments: &[ContourSegment],
    tolerance: f64,
) -> Result<Vec<ContourChain>, String> {
    stitch_impl(segments, tolerance, false)
}

/// Preserve every nondegenerate segment, stopping chains at junctions with
/// more than two neighbors. Never chooses a pairing across ambiguous vertices.
/// Callers may render resulting open chains as independent geometry.
pub fn stitch_segments_split_junctions(
    segments: &[ContourSegment],
    tolerance: f64,
) -> Result<Vec<ContourChain>, String> {
    stitch_impl(segments, tolerance, true)
}

fn stitch_impl(
    segments: &[ContourSegment],
    tolerance: f64,
    split_junctions: bool,
) -> Result<Vec<ContourChain>, String> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err("tolerance must be finite and positive".into());
    }
    type K = (i64,i64);
    let mut by_level: BTreeMap<u64, Vec<(K,K)>> = BTreeMap::new();
    let mut coordinates: BTreeMap<(u64,K), [f64;2]> = BTreeMap::new();
    for seg in segments {
        if !seg.elevation.is_finite() { return Err("nonfinite elevation".into()); }
        let a = key(seg.start,tolerance).ok_or("invalid start point")?;
        let b = key(seg.end,tolerance).ok_or("invalid end point")?;
        if a == b { continue; }
        let lev = seg.elevation.to_bits();
        coordinates.entry((lev,a)).or_insert(seg.start);
        coordinates.entry((lev,b)).or_insert(seg.end);
        by_level.entry(lev).or_default().push((a,b));
    }
    let mut output = Vec::new();
    for (lev, edges) in by_level {
        let mut adj: BTreeMap<K, Vec<(usize,K)>> = BTreeMap::new();
        let mut unique = BTreeSet::new();
        let mut kept = Vec::new();
        for (a,b) in edges {
            let canon = if a < b {(a,b)} else {(b,a)};
            if !unique.insert(canon) { continue; }
            let id = kept.len();
            kept.push((a,b));
            adj.entry(a).or_default().push((id,b));
            adj.entry(b).or_default().push((id,a));
        }
        if let Some((node, neighbors)) = adj.iter().find(|(_, neighbors)| neighbors.len() > 2).filter(|_| !split_junctions) {
            let xy = coordinates[&(lev, *node)];
            return Err(format!(
                "ambiguous contour junction (degree {}) at elevation {:.6}, X {:.6}, Y {:.6}; check exact-level TIN vertices or snapping tolerance",
                neighbors.len(), f64::from_bits(lev), xy[0], xy[1]
            ));
        }
        let mut visited = vec![false; kept.len()];
        // Open chains first, then cycles. Sorted map makes output deterministic.
        // A junction can have several independent outgoing chains. Revisit
        // each non-degree-two node for every incident edge before beginning
        // isolated degree-two cycles. Starting cycles too early could leave
        // two artificial open ends at an ordinary degree-two vertex.
        let starts: Vec<K> = adj.iter().filter(|(_,n)|n.len()!=2)
            .flat_map(|(k,n)|std::iter::repeat(*k).take(n.len()))
            .chain(adj.iter().filter(|(_,n)|n.len()==2).map(|(k,_)|*k))
            .collect();
        for start in starts {
            if !adj[&start].iter().any(|(id,_)| !visited[*id]) { continue; }
            let mut path = vec![start];
            let mut curr = start;
            loop {
                let next = adj[&curr].iter().filter(|(id,_)|!visited[*id])
                    .min_by_key(|(_,other)|*other).copied();
                let Some((id,other)) = next else {break};
                visited[id]=true;
                path.push(other);
                curr=other;
                if curr==start || (split_junctions && adj[&curr].len() != 2) {break}
            }
            let closed = path.len()>2 && path.first()==path.last();
            output.push(ContourChain {
                elevation: f64::from_bits(lev),
                vertices: path.iter().map(|k|coordinates[&(lev,*k)]).collect(),
                closed,
            });
        }
    }
    Ok(output)
}

/// Merge exactly four open chains ending at a verified saddle, using a
/// caller-selected non-crossing pairing. Default contour generation never
/// calls this function. Preserves all original points and elevation.
pub fn pair_chains_at_saddle(
    chains: &[ContourChain],
    center: [f64; 2],
    level: f64,
    tolerance: f64,
    choice: char,
) -> Result<Vec<ContourChain>, String> {
    if !matches!(choice,'A'|'B') || !level.is_finite()
        || !tolerance.is_finite() || tolerance<=0.0
        || center.iter().any(|v|!v.is_finite()) {
        return Err("invalid saddle pairing settings".into());
    }
    let close=|a:[f64;2]| (a[0]-center[0]).hypot(a[1]-center[1])<=tolerance;
    let mut tips:Vec<(f64,usize,Vec<[f64;2]>)>=Vec::new();
    for (i,c) in chains.iter().enumerate() {
        if c.elevation!=level || c.closed || c.vertices.len()<2 {continue;}
        let first=close(c.vertices[0]);
        let last=close(*c.vertices.last().unwrap());
        if first==last { continue; } // both or neither: unsafe to pair
        let oriented:Vec<[f64;2]>=if last {c.vertices.clone()}
            else {c.vertices.iter().rev().copied().collect()};
        let p=oriented[oriented.len()-2];
        let dx=p[0]-center[0];let dy=p[1]-center[1];
        if !dx.is_finite() || !dy.is_finite() || dx.hypot(dy)<=tolerance {
            continue;
        }
        tips.push((dy.atan2(dx),i,oriented));
    }
    // A contour loop based at the saddle contributes two incident half-edges,
    // although it is one closed chain. The former four-open-tip rule rejected
    // this valid graph. Keep all original segments; never duplicate a link.
    if tips.len()==2 {
        let loops:Vec<(usize,&ContourChain)>=chains.iter().enumerate().filter(|(_,c)|
            c.elevation==level && c.closed && c.vertices.len()>=4 &&
            close(c.vertices[0]) && close(*c.vertices.last().unwrap())
        ).collect();
        if loops.len()==1 {
            let (loop_id,cycle)=loops[0];
            let loop_segments=cycle.vertices.len()-1;
            let a=&tips[0].2;
            let b=&tips[1].2;
            let mut joined=a.clone(); // free -> saddle
            if choice=='B' {
                joined.extend(cycle.vertices.iter().skip(1).copied()); // saddle -> loop -> saddle
            }
            joined.extend(b[..b.len()-1].iter().rev().copied()); // saddle -> free
            let chosen:BTreeSet<usize>=[tips[0].1,tips[1].1,loop_id].into_iter().collect();
            let mut result:Vec<ContourChain>=chains.iter().enumerate()
                .filter(|(i,_)|!chosen.contains(i)).map(|(_,c)|c.clone()).collect();
            result.push(ContourChain{elevation:level,vertices:joined,closed:false});
            if choice=='A' {result.push(cycle.clone());}
            let input_edges:usize=chains.iter().map(|c|c.vertices.len().saturating_sub(1)).sum();
            let output_edges:usize=result.iter().map(|c|c.vertices.len().saturating_sub(1)).sum();
            if input_edges!=output_edges || loop_segments<3 {
                return Err("saddle loop edge conservation failed".into());
            }
            return Ok(result);
        }
    }
    if tips.len()!=4 {
        let mut closed_at_center=0usize;
        let mut center_to_center=0usize;
        let mut interior_only=0usize;
        for c in chains.iter().filter(|c|c.elevation==level && c.vertices.len()>=2) {
            let a=close(c.vertices[0]);
            let b=close(*c.vertices.last().unwrap());
            if a && b {
                if c.closed {closed_at_center+=1;} else {center_to_center+=1;}
            } else if !a && !b && c.vertices[1..c.vertices.len()-1].iter().any(|&p|close(p)) {
                interior_only+=1;
            }
        }
        return Err(format!(
            "expected four eligible saddle branches, found {}; {} saddle-to-saddle chain(s), {} closed chain(s) touching saddle at ends, {} chain(s) with saddle in interior. No changes applied.",
            tips.len(),center_to_center,closed_at_center,interior_only
        ));
    }
    tips.sort_by(|a,b|a.0.total_cmp(&b.0));
    let indices=if choice=='A' {[(0,1),(2,3)]}
        else {[(1,2),(3,0)]};
    let mut merged=Vec::new();
    for (a,b) in indices {
        let mut points=tips[a].2.clone(); // free endpoint -> saddle
        let right=&tips[b].2;
        points.extend(right[..right.len()-1].iter().rev().copied());
        merged.push(ContourChain{elevation:level,vertices:points,closed:false});
    }
    let selected:BTreeSet<usize>=tips.iter().map(|t|t.1).collect();
    let mut output:Vec<ContourChain>=chains.iter().enumerate()
        .filter(|(i,_)|!selected.contains(i)).map(|(_,c)|c.clone()).collect();
    output.extend(merged);
    Ok(output)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn segment(a:[f64;2],b:[f64;2])->ContourSegment{
        ContourSegment{elevation:5.,start:a,end:b}
    }
    #[test]
    fn degree_four_junction_exhausts_all_arms_before_cycles(){
        let segments=vec![
            segment([0.,0.],[1.,0.]),segment([1.,0.],[2.,0.]),
            segment([0.,0.],[-1.,0.]),segment([-1.,0.],[-2.,0.]),
            segment([0.,0.],[0.,1.]),segment([0.,1.],[0.,2.]),
            segment([0.,0.],[0.,-1.]),segment([0.,-1.],[0.,-2.]),
        ];
        let chains=stitch_segments_split_junctions(&segments,1e-8).unwrap();
        assert_eq!(chains.len(),4);
        assert_eq!(chains.iter().map(|c|c.vertices.len()-1).sum::<usize>(),8);
        let center_ends=chains.iter().flat_map(|c|
            [c.vertices.first(),c.vertices.last()]
        ).flatten().filter(|p|**p==[0.,0.]).count();
        assert_eq!(center_ends,4);
    }
    #[test]
    fn joins_reversed_segments_into_open_chain(){
        let x=stitch_segments(&[
            segment([1.,0.],[0.,0.]),
            segment([2.,0.],[1.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert_eq!(x[0].vertices.len(),3);
        assert!(!x[0].closed);
    }
    #[test]
    fn identifies_closed_loop(){
        let x=stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([1.,0.],[1.,1.]),
            segment([1.,1.],[0.,1.]),
            segment([0.,1.],[0.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert!(x[0].closed);
        assert_eq!(x[0].vertices.len(),5);
    }
    #[test]
    fn ignores_reversed_duplicate_edge(){
        let x=stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([1.,0.],[0.,0.]),
        ],1e-8).unwrap();
        assert_eq!(x.len(),1);
        assert_eq!(x[0].vertices.len(),2);
    }
    #[test]
    fn refuses_branch_points(){
        assert!(stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([0.,0.],[0.,1.]),
            segment([0.,0.],[-1.,0.])
        ],1e-8).is_err());
    }
    #[test]
    fn branch_diagnostic_identifies_level_and_position(){
        let err = stitch_segments(&[
            segment([0.,0.],[1.,0.]),
            segment([0.,0.],[0.,1.]),
            segment([0.,0.],[-1.,0.])
        ],1e-8).unwrap_err();
        assert!(err.contains("degree 3"));
        assert!(err.contains("elevation 5.000000"));
        assert!(err.contains("X 0.000000"));
    }
    #[test]
    fn splits_degree_four_without_connecting_unrelated_branches() {
        let segs = [
            segment([-2.,0.],[-1.,0.]),
            segment([-1.,0.],[0.,0.]),
            segment([0.,0.],[1.,0.]),
            segment([0.,0.],[0.,1.]),
            segment([0.,0.],[0.,-1.]),
        ];
        let chains = stitch_segments_split_junctions(&segs,1e-8).unwrap();
        assert_eq!(chains.len(),4);
        assert_eq!(chains.iter().map(|c|c.vertices.len()-1).sum::<usize>(),5);
        assert!(chains.iter().all(|c|!c.closed));
    }
    #[test]
    fn selected_saddle_pairing_preserves_segment_count() {
        let c=[0.,0.];
        let branches=vec![
            ContourChain{elevation:5.,vertices:vec![[1.,0.],c],closed:false},
            ContourChain{elevation:5.,vertices:vec![c,[0.,1.]],closed:false},
            ContourChain{elevation:5.,vertices:vec![[-1.,0.],c],closed:false},
            ContourChain{elevation:5.,vertices:vec![c,[0.,-1.]],closed:false}
        ];
        for choice in ['A','B'] {
            let combined=pair_chains_at_saddle(&branches,c,5.,1e-8,choice).unwrap();
            assert_eq!(combined.len(),2);
            assert_eq!(combined.iter().map(|p|p.vertices.len()-1).sum::<usize>(),4);
            assert!(combined.iter().all(|p|p.elevation==5.));
        }
    }
    #[test]
    fn pairing_rejects_insufficient_branches_without_changes() {
        let branch=ContourChain{elevation:5.,vertices:vec![[1.,0.],[0.,0.]],closed:false};
        assert!(pair_chains_at_saddle(&[branch.clone()],[0.,0.],5.,1e-8,'A').is_err());
        assert_eq!(branch.vertices.len(),2);
    }
    #[test]
    fn saddle_loop_reports_missing_free_branches_without_mutating_input() {
        let c=[0.,0.];
        let chains=vec![
            ContourChain{elevation:5.,vertices:vec![[1.,0.],c],closed:false},
            ContourChain{elevation:5.,vertices:vec![[0.,1.],c],closed:false},
            ContourChain{elevation:5.,vertices:vec![c,[-1.,0.],[0.,-1.],c],closed:false},
        ];
        let saved=chains.clone();
        let error=pair_chains_at_saddle(&chains,c,5.,1e-8,'A').unwrap_err();
        assert!(error.contains("found 2"));
        assert!(error.contains("1 saddle-to-saddle chain"));
        assert_eq!(chains,saved);
    }
    #[test]
    fn closed_saddle_loop_and_two_open_arms_preserve_all_edges() {
        let c=[0.,0.];
        let chains=vec![
            ContourChain{elevation:5.,vertices:vec![[2.,0.],[1.,0.],c],closed:false},
            ContourChain{elevation:5.,vertices:vec![c,[0.,1.],[0.,2.]],closed:false},
            ContourChain{elevation:5.,vertices:vec![c,[-1.,0.],[-1.,-1.],[0.,-1.],c],closed:true}
        ];
        let saved=chains.clone();
        for choice in ['A','B'] {
            let out=pair_chains_at_saddle(&chains,c,5.,1e-8,choice).unwrap();
            assert_eq!(out.iter().map(|x|x.vertices.len()-1).sum::<usize>(),8);
            assert!(out.iter().all(|x|x.elevation==5.));
            assert_eq!(out.len(),if choice=='A'{2}else{1});
            assert_eq!(chains,saved);
        }
    }
    #[test]
    fn rejects_invalid_tolerance(){
        assert!(stitch_segments(&[],0.).is_err());
    }
}
