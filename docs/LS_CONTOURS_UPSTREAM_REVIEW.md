# LS_CONTOURS upstream contribution plan (Kevin's review, 2026-10-22)

Status: preparation only. **Do not open a PR yet**: this branch is still based on the fork's 0.4.19 main, not the upstream 0.4.24/v2026.42 tree. The experimental development branches with profiles and grades must remain independent.

## Fixed scope
* Engine-only contour geometry in `crates/landsurvey/`, using existing `Surface`, triangle indexing and `surface.rs` isoline/cut-fill helpers. No CAD dependencies in engine.
* Stitch unordered segments into deterministic open/closed chains; test saddle, exact-level vertices, flat edges, degenerates, topology and length.
* CAD: `Polyline3D`, layers `LS-CONTOUR-MAJ`, `LS-CONTOUR-MIN`, labels `LS-CONTOUR-LABEL`.
* XDATA: reuse `LANDSURVEY_SURFACE` kind `CONTOUR`, optionally elevation as a third value, preserving the name+kind prefix.
* CAD insertions via `src/batch.rs::Batch`, not per-entity round trips; remain under the 30-second host dispatch deadline.
* Resolve surfaces via `resolve_named_surface`; add Surface ribbon action with a command arm.
* Tests: engine unit cases, independent numeric golden, CAD XDATA/DWG save-reopen and regeneration tests.
* Major/minor command UX: `LS_CONTOURS <surface-name> [minor-interval=0.5] [major-every=5]`.
* No smoothing, profiles, grades, or other exploratory features in the initial upstream PR.

## Gaps in this existing experimental branch
* `contours.rs` currently has a distinct half-open edge-level slicer; reconcile with `surface.rs::tri_cut_fill_seg` and preserve its exact-vertex behavior.
* Focused branch is based on host v2026.38/0.4.19. It must be rebased/ported to upstream v2026.42/0.4.24 before any PR.
* Experimental host drawing uses `LS-CONTOUR-MAJOR/MINOR` and per-entity `tag_surface`; adjust to upstream layers and `Batch`.
* Independent golden counts and contour lengths are not yet externally verified.

## Merge gates
1. Diff strictly limited to contours-related changes against upstream main.
2. `cargo test --workspace --locked` and `cargo build --release --locked` pass on v2026.42/Rust 1.99.1.
3. Geometry, XDATA, DWG round-trip and no-duplicate regeneration tests pass.
4. Share La Piedra demo as illustrative fixture only, not as independently certified golden.
