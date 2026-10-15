//! Export generated TIN contours as R12 DXF 3D polylines.
//! No host API usage: safe for headless CLI and tests.
use crate::contour_levels::ClassifiedContour;
use crate::dxf::DxfBuilder;

/// Render major and minor contour chains in separate layers.
/// Returns valid R12 ASCII DXF; Z comes from each contour elevation.
pub fn contours_to_dxf(contours: &[ClassifiedContour]) -> Result<String, String> {
    let mut dxf = DxfBuilder::new();
    dxf.add_layer("LS-CONTOUR-MAJOR", 2);
    dxf.add_layer("LS-CONTOUR-MINOR", 8);
    for contour in contours {
        let chain = &contour.chain;
        if !chain.elevation.is_finite() || chain.vertices.len() < 2 ||
            chain.vertices.iter().flatten().any(|v| !v.is_finite()) {
            return Err("invalid contour geometry".into());
        }
        let layer = if contour.major {"LS-CONTOUR-MAJOR"} else {"LS-CONTOUR-MINOR"};
        let points: Vec<[f64;3]> = chain.vertices.iter()
            .map(|p| [p[0],p[1],chain.elevation]).collect();
        dxf.polyline3d(layer, &points, chain.closed);
    }
    Ok(dxf.build())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::surface::Surface;
    use crate::contour_levels::generate_contours;
    #[test]
    fn exports_half_meter_contours_as_true_3d_polylines() {
        let surface=Surface{nodes:vec![[355100.,928700.,0.],
            [355110.,928700.,5.],[355100.,928710.,0.]],triangles:vec![[0,1,2]]};
        let result=generate_contours(&surface,0.5,5,0.,1e-6).unwrap();
        let dxf=contours_to_dxf(&result).unwrap();
        assert_eq!(dxf.matches("0\nPOLYLINE\n").count(),9);
        assert!(dxf.contains("8\nLS-CONTOUR-MAJOR\n"));
        assert!(dxf.contains("8\nLS-CONTOUR-MINOR\n"));
        assert!(dxf.contains("30\n2.500000\n"));
        assert!(dxf.contains("10\n3551"));
        assert_eq!(dxf.matches("0\nSEQEND\n").count(),9);
    }
    #[test]
    fn invalid_geometry_is_rejected(){
        use crate::contour_chains::ContourChain;
        let bad=ClassifiedContour{major:true,chain:ContourChain{
            elevation:f64::NAN,vertices:vec![[0.,0.],[1.,1.]],closed:false}};
        assert!(contours_to_dxf(&[bad]).is_err());
    }
}
