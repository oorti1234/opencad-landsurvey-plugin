# Focused CAD integration plan (not yet implemented)

The headless contour engine and its DXF regressions are already isolated on this
branch. Do **not** copy the complete `src/dispatch.rs` from
`release/contours-integration-candidate`: that file includes unrelated
LandXML, field-to-finish, UI, and host compatibility work.

## Source inventory

Validated implementation reference: `release/contours-integration-candidate`.

- `fn contours_command` in `src/dispatch.rs`: roughly 142 lines.
- `fn contour_quality_command`: roughly 84 lines.
- `fn contour_label_position`: roughly 23 lines.
- `fn resolve_named_surface` and `fn surface_record` already exist in
  the base `main` and should be reused where compatible.
- `src/ribbon.rs` and command routing require small, selective patches.

## Safe integration order

1. Compare the host/codec ABI and `Batch` implementation against `main`.
   The release branch imports `codec::{Polyline3D, Color}` whereas the
   unmodified base imports `acadrust` types. Verify compatibility before
   copying CAD helpers.
2. Add only the `LS_CONTOURS` route and associated CAD entity generation,
   XDATA kind, replacement, label positioning, and undo transaction.
3. Add only the read-only `LS_CONTOUR_QA` route. Keep saddle diagnostics and
   smoothing UI optional follow-up commits, not a prerequisite.
4. Add the minimal ribbon entry and usage documentation, avoiding unrelated
   button/menu changes.
5. Run `cargo test --workspace --locked` and
   `cargo build --release --locked` on this exact branch.
6. Install only after build success, and test on copies of two DWGs.
7. Re-run contours twice; check no duplicates, accurate XDATA lookup, undo/redo,
   DWG save/reopen, unchanged TIN, and the two QA baselines documented under
   `docs/CONTOURS_INTEGRATION_CANDIDATE.md` on the release branch.

## Observed QA baselines (full integration branch)

- La Piedra: 41 chains, 2555 segments, 0 unexplained interior endpoints.
- Second PNEZD: 23 chains, 531 segments, 0 unexplained interior endpoints.

Do not interpret the above as tests of this focused branch's CAD integration:
the focused branch currently contains the headless engine only.
