# Contributing to wall-ctf

Thank you for your interest in contributing to wall-ctf!

## How to report issues

If you find a bug or have a feature request, please open an issue on GitHub:
https://github.com/valeriolobrano/wall-ctf/issues

When reporting a bug, please include:
- A minimal wall definition (JSON) that reproduces the problem
- The parameters used (n_roots, n_coefficients, sampling_time)
- The error message or unexpected output
- Your Python and wall-ctf versions (`python --version` and `pip show wall-ctf`)

## How to contribute code

1. Fork the repository
2. Create a feature branch (`git checkout -b feature/my-improvement`)
3. Make your changes
4. Run the tests (`uv run pytest`)
5. Commit and push
6. Open a Pull Request

Please ensure that:
- All existing tests pass
- New features include appropriate tests
- Code follows the existing style (no additional linting rules imposed)

## How to seek support

For questions about usage, open a GitHub Discussion or issue.
For academic collaboration, contact: valerio.lobrano@unipa.it

## Citation

If you use wall-ctf in your research, please cite the publications listed
in the [README](README.md#how-to-cite).
