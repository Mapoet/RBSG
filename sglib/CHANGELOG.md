# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.1.0] - 2025-05-19

### Added
- Initial release of the Robust Savitzky-Golay Filter Library
- Basic SG filter implementation
- Coefficient calculation module
- Standard SG filtering implementation
- Robust filtering with configurable outlier detection
- Command-line test program
- Python visualization tools
- Documentation and examples

## [0.2.0] - 2025-05-20

### Added
- Added configurable rate parameter for outlier detection sensitivity
- Command-line option `-r/--rate` to adjust outlier detection threshold
- Updated documentation to explain the rate parameter usage

### Changed
- Improved outlier detection algorithm with iterative refinement
- Enhanced robustness of window border handling
- Updated examples to demonstrate rate parameter usage

### Fixed
- Corrected coefficient normalization in edge cases
- Fixed memory allocation issues in large datasets 