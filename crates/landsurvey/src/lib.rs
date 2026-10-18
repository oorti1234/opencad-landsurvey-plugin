//! Headless COGO + PNEZD survey geometry for Open CAD Studio.
//!
//! Layer-C engine (see `docs/plugin-architecture.md`): depends only on `std`
//! and `serde` — no `iced`, `acadrust`, or host types — so it builds unchanged
//! for CLI/WASM. Its in-app consumer is the `opencad-landsurvey-plugin` cdylib,
//! which wraps these functions in ribbon/command/XDATA glue.

pub mod cogo;
pub mod contours;
pub mod contour_levels;
pub mod contour_topology;
pub mod contour_quality;
pub mod contour_smoothing;
pub mod contour_smoothing_qa;
pub mod saddle_pairing;
pub mod contour_dxf;
pub mod contour_chains;
pub mod tin_spatial_index;
pub mod dxf;
pub mod landxml;
pub mod plan;
pub mod pnezd;
pub mod resection;
pub mod surface;
pub mod transform;
pub mod viz;
