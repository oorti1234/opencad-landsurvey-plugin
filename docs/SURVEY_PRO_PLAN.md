# OpenCAD Survey Pro — compatibility and development plan

Target environment: Ubuntu 26.04.1 LTS (x86_64), Open CAD Studio v2026.40.1 AppImage.

## Baseline from the existing fork
- Land Survey plugin v0.4.19.
- `host-build.json` pins Open CAD Studio **v2026.38** (commit `0d023d267bc5b7afeca3b54e98875b0efd4f3926`).
- `rust-toolchain.toml` pins Rust 1.98.1.
- `plugin.toml` declares `api_version = 5`.
- Current functionality documented in `PLUGIN.md`: `LS_PNEZD`, `LS_SURFACE`, COGO, transforms, volumes.

**Do not claim this plugin is compatible with v2026.40.1 until build and smoke tests confirm it.**

## Stage 1: establish compatibility
1. Inspect the source and lockfile of the Open CAD Studio v2026.40.1 release.
2. Compare API_VERSION, ocs_plugin_api and acadrust/cadcodec revisions, compiler and ABI requirements.
3. Use the project's existing repin scripts and gates to update all relevant pins consistently.
4. Build with `cargo build --release --locked`; run `cargo test --workspace --locked`.
5. Build Linux x86_64 release assets in CI and confirm the manifest metadata.
6. Smoke-test installation through Plugin Manager in Open CAD Studio v2026.40.1. Confirm Land Survey ribbon and `LS_HELLO`, `LS_PNEZD`, `LS_SURFACE` on a copy of test drawings.
7. Never replace a working plugin installation until the new release passes all tests.

## Stage 2: contours
Add a host-independent contour generator to `crates/landsurvey`. Interpolate isolines from TIN triangles, join segments into lines with controlled tolerances, handle vertex-on-level and flat-edge cases, and provide deterministic output. Validate on planar slopes, saddles, duplicates and degenerate triangles. Expose optional 0.50 m minor and 2.50 m major contours and DXF export after the core algorithm passes tests.

## Preservation rules
- All changes occur on branch `desarrollo-survey-pro`, not `main`.
- Keep the existing commands and data formats intact.
- Keep an untouched copy of the v2026.36 and v2026.40.1 AppImages.
- Treat survey data and terrain examples as synthetic, not surveyed field observations.
- Record test outputs before marking work as released.
