# Changelog - MICE CloWM Workflow

## Unreleased (Work in Progress)

## 0.1.4 - 2026-07-30

### Added
- First draft of `clowm/README.md`, `clowm/usage.md` and `clowm/output.md` 

## 0.1.3 - 2026-07-30

### Added
- S3 URI validation for output directory path
- Cleaning up artefact '.' directory created as part of workflow execution

### Changed
- Workflow now modularised to support RO-Crate generation

## 0.1.2 - 2026-07-24

### Fixed

- Handling of output directory; no longer wipes working dir if supplying "./" to --out_dir

## 0.1.1 - 2026-06-10

### Added

- First draft of parameter schema
- Docker container now in GHCR; updated `main.nf` accordingly

## 0.1.0 - 2026-06-05

### Added

- Initial Docker container with `mice`, `ggcat` and `gfa2gff`
- Initial version of NF workflow, providing just `mice`
