# Contributing to SG-Lib

Thank you for considering contributing to the Robust Savitzky-Golay Filter Library (SG-Lib). This document provides guidelines and instructions for contributing to this project.

## Code of Conduct

Please be respectful and considerate of others. We aim to foster an inclusive and welcoming community.

## How to Contribute

### Reporting Bugs

If you find a bug, please create an issue with the following information:

- A clear, descriptive title
- Steps to reproduce the bug
- Expected behavior
- Actual behavior
- Any additional context (OS, compiler version, etc.)

### Suggesting Enhancements

For feature requests or enhancements:

- Clearly describe the feature/enhancement
- Explain why it would be useful
- Provide examples of how it would be used

### Pull Requests

1. Fork the repository
2. Create a new branch from `main`
3. Make your changes
4. Run tests to ensure your changes don't break existing functionality
5. Submit a pull request

### Coding Standards

#### Fortran Code

- Follow modern Fortran standards (Fortran 2008+)
- Use explicit typing and `implicit none`
- Indent with 2 spaces
- Use lowercase for variables and procedures
- Follow the existing code style and structure
- Add appropriate comments and documentation

#### Python Scripts

- Follow PEP 8 style guide
- Use descriptive variable and function names
- Add docstrings for functions and modules

## Development Setup

1. Clone the repository
2. Install dependencies:
   - Fortran compiler (gfortran 8+ recommended)
   - OpenBLAS and LAPACK libraries
   - Python 3.6+ for scripts (optional)

3. Build the library:
   ```bash
   make
   ```

4. Run tests:
   ```bash
   make test
   ```

## Release Process

- Versioning follows [Semantic Versioning](https://semver.org/)
- Changes for each release are documented in CHANGELOG.md

## License

By contributing to this project, you agree that your contributions will be licensed under the project's [MIT License](LICENSE).

## Questions?

If you have any questions about contributing, please open an issue and we'll be happy to help! 