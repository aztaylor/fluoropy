# Contributing

## Setup

```bash
pip install -e ".[dev,viz]"
pre-commit install
```

`pre-commit install` runs black, isort and flake8 on every commit, at the
same pinned versions CI uses.

## Before pushing

```bash
pytest
black --check fluoropy tests
isort --check-only fluoropy tests
flake8 fluoropy tests
```

CI runs all four. To fix formatting rather than just check it, drop
`--check` / `--check-only`.

## Conventions

See [ARCHITECTURE.md](ARCHITECTURE.md) for where code goes and the
conventions it follows: array axes, roles, well IDs, logging, and failing
loudly instead of returning plausible wrong numbers.

Bulk reformatting commits go in `.git-blame-ignore-revs` so they stay out of
`git blame`.
