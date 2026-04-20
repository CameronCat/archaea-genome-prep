# Contributing to archaea-genome-prep

Thank you for your interest in contributing!

## Reporting bugs

Please open an issue on GitHub with:
- A description of the bug
- A minimal reproducible example (FASTA file or synthetic sequence)
- Your Python version and operating system

## Suggesting features

Open an issue describing:
- The feature you would like
- The biological use case it addresses
- Any relevant tools or literature

## Contributing code

1. Fork the repository
2. Create a new branch: `git checkout -b my-feature`
3. Make your changes
4. Run the tests: `pytest`
5. Push your branch and open a pull request

## Development setup

```bash
git clone https://github.com/CameronPiepkorn/archaea-genome-prep
cd archaea-genome-prep
pip install -e ".[dev]"
pytest
```

## Code style

- Follow PEP 8
- Add docstrings to new functions
- Add tests for new functionality
- No hard dependencies — keep the library pure stdlib only
- Prokka integration must remain a soft dependency (detected at runtime)

## Questions

Open an issue on GitHub.
