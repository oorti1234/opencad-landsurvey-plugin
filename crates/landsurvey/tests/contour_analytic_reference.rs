//! Analytic contour reference independent of the production slicing/stitching code.
//!
//! On the rectangle 0 <= x,y <= 2, the planar TIN has elevation z=x.
//! Therefore the exact isoline z=h is x=h, 0 <= y <= 2, for 0<h<2.
//! Every such level has precisely one open contour of length 2.
//! These expectations are derived algebraically, never from the implementation.
//!
//! This is an independent analytic reference, not an external Civil 3D golden.
use landsurvey::contour_chains::stitch_segments;
use landsurvey::contours::segments_at_level;
use landsurvey::surface::Surface;

#[test]
fn planar_multilevel_contours_match_exact_reference() {
    let surface = Surface {
        nodes: vec![
            [0., 0., 0.], [2., 0., 2.],
            [0., 2., 0.], [2., 2., 2.],
        ],
        triangles: vec![[0, 1, 2], [1, 3, 2]],
    };

    for level in [0.5, 1.0, 1.5] {
        let segments = segments_at_level(&surface, level);
        let chains = stitch_segments(&segments, 1e-9)
            .expect("an analytic plane cannot have ambiguous branches");

        assert_eq!(chains.len(), 1, "one connected contour at z={level}");
        let chain = &chains[0];
        assert!(!chain.closed, "isoline must reach opposite boundary edges");
        assert_eq!(chain.elevation, level);
        assert!(chain.points.len() >= 2);
        assert!(chain.points.iter().all(|p| (p[0] - level).abs() < 1e-10),
            "on z=x, every contour vertex must have x equal to its level");

        let first = chain.points.first().unwrap();
        let last = chain.points.last().unwrap();
        assert!((first[1].min(last[1]) - 0.0).abs() < 1e-10);
        assert!((first[1].max(last[1]) - 2.0).abs() < 1e-10);

        let length: f64 = chain.points.windows(2)
            .map(|ends| (ends[1][0] - ends[0][0])
                .hypot(ends[1][1] - ends[0][1]))
            .sum();
        assert!((length - 2.0).abs() < 1e-10,
            "exact analytic contour length is 2 at z={level}, got {length}");
    }
}
