# LS_CONTOURS — integration candidate validation

This branch consolidates the experimental TIN contour generator, CAD 3D polylines,
major labels, safe regeneration, QA, saddle diagnostics and protected smoothing.

## Commands

- `LS_SURFACE /absolute/path/to/file.txt` — import and retain a named TIN, drawing its edges.
- `LS_CONTOURS SURFACE_NAME 0.5 5` — contour minor interval 0.5, majors every 5 levels.
- `LS_CONTOUR_QA SURFACE_NAME 0.5 5` — inspect boundary, junction and unexplained endpoints.
- `LS_CONTOUR_SMOOTH SURFACE_NAME 0.5 1 0.50 LOCAL 0.010` — independent preview in a separate layer, sampled XY and Z limits.
- `LS_SADDLE_ASSEMBLE SURFACE_NAME LEVEL A` — experimental choice of saddle connectivity; not a definitive topographic interpretation.
- `LS_SADDLE_CLEAR SURFACE_NAME` — remove tagged experimental saddle entities.

## Reproduced UI observations

| QA metric | La Piedra synthetic terrain | Second PNEZD surface (grid-priority) |
|---|---:|---:|
| contour chains | 41 | 23 |
| open / closed | 25 / 16 | 18 / 5 |
| contour edges | 2555 | 531 |
| boundary endpoints | 48 | 34 |
| ambiguous endpoints | 2 | 2 |
| unexplained interior endpoints | 0 | 0 |
| zero-length / duplicate edges | 0 / 0 | 0 / 0 |

The second surface previously produced 25 chains and four unexplained interior
endpoints. The adjusted chain traversal produced 23 chains and zero unexplained
endpoints while retaining all 531 edges. One exact-level saddle remains and is
not resolved automatically.

## Release gate (run locally)

```bash
cargo test --workspace --locked
cargo build --release --locked
python3 .github/scripts/build_metadata.py --output dist/plugin.toml
```

Test on **copies** of both DWGs. Confirm `LS_CONTOUR_QA` on each surface,
repeat `LS_CONTOURS` and confirm safe replacement rather than duplication,
save/reopen DWG, test CAD undo and redo, and compare original contours with the
smoothing preview. Verify host/plugin ABI compatibility.

Known limitations: the smoothing residuals and crossing checks are sampled or
tolerance-based and are not topographic certification; saddle choice remains
experimental; point-duplicate source data require surveyor review. This candidate
has not yet been accepted upstream.
