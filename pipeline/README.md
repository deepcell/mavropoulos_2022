# askCell v0.0.1

This package provides software to analyze molecular data based on multiple
types of sequencing data.

## Developer installation

First, set up your base Python virtual environment however you wish (using conda,
pyenv, virtualenv, etc.).

Install development pre-requirements:

```bash
$ pip install -U pip poetry poetry-dynamic-versioning pre-commit
$ poetry config virtualenvs.create false  # if you wish to manage your own virtualenvs
```

Clone the repository, set up pre-commit hooks, and install from your
preferred project directory:

```bash
$ cd pipeline
$ pre-commit install
$ poetry install
```

Run tests to ensure everything is working:

```base
$ pytest
```

## Documentation

To generate the user and developer manuals for this package, go to this
package's base directory and run:

    make html

The manuals will be built and placed in build/html. The entry point is
build/html/index.html
