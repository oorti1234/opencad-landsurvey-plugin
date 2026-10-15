//! DXF contour regression tests: topology, coordinate preservation and layers.
use landsurvey::contour_chains::ContourChain;
use landsurvey::contour_dxf::contours_to_dxf;
use landsurvey::contour_levels::ClassifiedContour;

fn pairs(dxf: &str) -> Vec<(&str, &str)> {
    let lines: Vec<_> = dxf.lines().collect();
    assert_eq!(lines.len() % 2, 0, "DXF group code/value mismatch");
    lines.chunks_exact(2).map(|p|(p[0],p[1])).collect()
}

#[test]
fn export_preserves_utm_coordinates_and_z_per_vertex() {
    let chain = ContourChain {
        elevation: 42.5,
        vertices: vec![[355189.547,928708.151],[355200.021,928719.4]],
        closed: false,
    };
    let dxf = contours_to_dxf(&[ClassifiedContour{chain,major:true}]).unwrap();
    let groups = pairs(&dxf);
    assert!(groups.contains(&("10","355189.547000")));
    assert!(groups.contains(&("20","928708.151000")));
    assert_eq!(groups.iter().filter(|&&(code,val)| code=="30" && val=="42.500000").count(), 2);
    assert!(groups.contains(&("8","LS-CONTOUR-MAJOR")));
}

#[test]
fn closed_chain_has_closed_flag_and_no_duplicate_final_vertex() {
    let chain = ContourChain {
        elevation: 40.0,
        vertices: vec![[0.,0.],[1.,0.],[1.,1.],[0.,0.]],
        closed: true,
    };
    let dxf = contours_to_dxf(&[ClassifiedContour{chain,major:false}]).unwrap();
    let groups = pairs(&dxf);
    let polylines: Vec<_> = groups.iter().enumerate()
       .filter(|(_,p)| **p==("0","POLYLINE")).map(|(i,_)|i).collect();
    assert_eq!(polylines.len(),1);
    assert!(groups[polylines[0]..].iter().take(5).any(|p| *p==("70","9")));
    assert_eq!(groups.iter().filter(|p| **p==("0","VERTEX")).count(),3);
    assert!(groups.contains(&("8","LS-CONTOUR-MINOR")));
}

#[test]
fn empty_export_has_valid_structure() {
    let dxf=contours_to_dxf(&[]).unwrap();
    assert!(dxf.starts_with("0\nSECTION\n"));
    assert!(dxf.ends_with("0\nEOF\n"));
    assert!(!dxf.contains("0\nPOLYLINE\n"));
    pairs(&dxf);
}
