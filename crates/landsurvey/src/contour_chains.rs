//! Join independent TIN intersection segments into continuous contour chains.
//!
//! Branching junctions are reported, not silently resolved: saddle pairing
//! requires an explicit topology-aware policy in a later contribution step.
use crate::contours::ContourSegment;
use std::collections::BTreeMap;

#[derive(Debug, Clone, PartialEq)]
pub struct ContourChain {
    pub elevation: f64,
    pub points: Vec<[f64; 2]>,
    pub closed: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub enum ChainError {
    InvalidTolerance,
    InvalidSegment,
    Branch { elevation: f64, degree: usize },
}

fn vertex_key(point: [f64; 2], tolerance: f64) -> Option<(i64, i64)> {
    let x = (point[0] / tolerance).round();
    let y = (point[1] / tolerance).round();
    if !x.is_finite()
        || !y.is_finite()
        || x < i64::MIN as f64
        || x >= i64::MAX as f64
        || y < i64::MIN as f64
        || y >= i64::MAX as f64
    {
        return None;
    }
    Some((x as i64, y as i64))
}

/// Join all edges at a single elevation, with endpoint tolerance in XY.
///
/// An endpoint may meet only two distinct segments. Higher-degree junctions
/// return an error rather than inventing an arbitrary or crossing pairing.
/// A `closed` chain has no repeated last point.
pub fn stitch_segments(
    segments: &[ContourSegment],
    tolerance: f64,
) -> Result<Vec<ContourChain>, ChainError> {
    if !tolerance.is_finite() || tolerance <= 0.0 {
        return Err(ChainError::InvalidTolerance);
    }
    if segments.is_empty() {
        return Ok(Vec::new());
    }
    let level = segments[0].elevation;
    let mut vertices = BTreeMap::<(i64, i64), usize>::new();
    let mut coords = Vec::<[f64; 2]>::new();
    let mut edges = Vec::<(usize, usize)>::new();
    for segment in segments {
        if !level.is_finite() || segment.elevation != level {
            return Err(ChainError::InvalidSegment);
        }
        let mut ends = [0usize; 2];
        for (i, point) in [segment.start, segment.end].into_iter().enumerate() {
            let key = vertex_key(point, tolerance).ok_or(ChainError::InvalidSegment)?;
            ends[i] = *vertices.entry(key).or_insert_with(|| {
                let id = coords.len();
                coords.push(point);
                id
            });
        }
        if ends[0] == ends[1] {
            return Err(ChainError::InvalidSegment);
        }
        edges.push((ends[0], ends[1]));
    }
    let mut neighbors = vec![Vec::<usize>::new(); coords.len()];
    for (index, &(a, b)) in edges.iter().enumerate() {
        neighbors[a].push(index);
        neighbors[b].push(index);
    }
    for incident in &neighbors {
        if incident.len() > 2 {
            return Err(ChainError::Branch {
                elevation: level,
                degree: incident.len(),
            });
        }
    }
    let mut visited = vec![false; edges.len()];
    let mut chains = Vec::new();
    let mut seeds: Vec<usize> = (0..coords.len())
        .filter(|&id| neighbors[id].len() == 1)
        .collect();
    seeds.extend((0..coords.len()).filter(|&id| neighbors[id].len() == 2));
    for start in seeds {
        if neighbors[start].iter().all(|&edge| visited[edge]) {
            continue;
        }
        let mut points = vec![coords[start]];
        let mut current = start;
        let mut closed = false;
        loop {
            let edge = neighbors[current].iter().copied().find(|&e| !visited[e]);
            let Some(edge) = edge else { break };
            visited[edge] = true;
            let (a, b) = edges[edge];
            current = if a == current { b } else { a };
            if current == start {
                closed = true;
                break;
            }
            points.push(coords[current]);
        }
        chains.push(ContourChain {
            elevation: level,
            points,
            closed,
        });
    }
    debug_assert!(visited.iter().all(|&v| v));
    Ok(chains)
}

#[cfg(test)]
mod tests {
    use super::*;
    fn segment(a: [f64; 2], b: [f64; 2]) -> ContourSegment {
        ContourSegment {
            elevation: 10.,
            start: a,
            end: b,
        }
    }

    #[test]
    fn reversed_edges_join_into_open_chain() {
        let edges = [segment([1., 0.], [0., 0.]), segment([1., 0.], [2., 0.])];
        let chains = stitch_segments(&edges, 1e-6).unwrap();
        assert_eq!(chains.len(), 1);
        assert!(!chains[0].closed);
        assert_eq!(chains[0].points.len(), 3);
    }

    #[test]
    fn closed_ring_does_not_repeat_first_vertex() {
        let edges = [
            segment([0., 0.], [1., 0.]),
            segment([1., 1.], [0., 1.]),
            segment([1., 0.], [1., 1.]),
            segment([0., 1.], [0., 0.]),
        ];
        let chains = stitch_segments(&edges, 1e-6).unwrap();
        assert_eq!(chains.len(), 1);
        assert!(chains[0].closed);
        assert_eq!(chains[0].points.len(), 4);
    }

    #[test]
    fn disconnected_segments_form_two_chains() {
        let edges = [segment([0., 0.], [1., 0.]), segment([5., 0.], [6., 0.])];
        assert_eq!(stitch_segments(&edges, 1e-6).unwrap().len(), 2);
    }

    #[test]
    fn branching_junction_is_rejected_explicitly() {
        let edges = [
            segment([0., 0.], [1., 0.]),
            segment([1., 0.], [2., 0.]),
            segment([1., 0.], [1., 1.]),
        ];
        assert_eq!(
            stitch_segments(&edges, 1e-6),
            Err(ChainError::Branch {
                elevation: 10.,
                degree: 3
            })
        );
    }

    #[test]
    fn four_arm_saddle_is_reported_instead_of_cross_connected() {
        let edges = [
            segment([-1., 0.], [0., 0.]),
            segment([0., 0.], [1., 0.]),
            segment([0., -1.], [0., 0.]),
            segment([0., 0.], [0., 1.]),
        ];
        assert_eq!(
            stitch_segments(&edges, 1e-8),
            Err(ChainError::Branch {
                elevation: 10.,
                degree: 4
            })
        );
    }

    #[test]
    fn mismatched_elevation_is_rejected() {
        let edges = [
            segment([0., 0.], [1., 0.]),
            ContourSegment {
                elevation: 11.,
                start: [1., 0.],
                end: [2., 0.],
            },
        ];
        assert_eq!(
            stitch_segments(&edges, 1e-8),
            Err(ChainError::InvalidSegment)
        );
    }

    #[test]
    fn invalid_input_is_rejected() {
        assert_eq!(stitch_segments(&[], 0.), Err(ChainError::InvalidTolerance));
        assert_eq!(
            stitch_segments(&[segment([0., 0.], [0., 0.])], 1e-6),
            Err(ChainError::InvalidSegment)
        );
    }
}
