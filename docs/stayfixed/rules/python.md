# Python

Recommended defaults for this repository's Python, not hard rules: where the project already
does something else on purpose, the project wins. Each line is here because agents get it wrong
without being told.

## Before the first command

- Run commands through the package manager whose lockfile is committed (`uv run`,
  `poetry run`, `pdm run`, `pipenv run`); do not switch tools or add a second lockfile.
- After editing dependencies, re-lock and commit the lockfile in the same change; CI refuses a
  stale lock (`uv sync --locked`, `poetry check --lock`, `pdm install --check`,
  `pipenv install --deploy`).
- Run the test command CI runs: a `-m` on the command line replaces the one in `addopts`, and a
  misspelt marker in `-m "not x"` deselects nothing and passes.

## Package manager and lockfile

- An application commits its lockfile. A library may run unlocked on purpose, and then lists
  the lockfile in `.gitignore`, because the package manager writes it again on the next run.
- `pylock.toml` is the standard lock format. A project locked with `pip lock` commits it as its
  lock; uv and PDM can export it beside their own lock, and still install from their own.

## Python version

- `requires-python` declares a minimum and no upper bound; tested versions go in classifiers.
- Where the package manager supports it, development dependencies go in `[dependency-groups]`
  rather than in an extra that installs with the package.

## Lint, format, types

- If the project uses Ruff, an explicit `select` replaces its defaults entirely; extend instead
  of replacing unless the narrower set is deliberate.
- One type checker, configured and pinned, whichever the project chose. Do not add a second
  one, and do not "fix" a finding by switching checkers.
- Declare enums a type checker must see with class syntax. The functional `Enum("Name", ...)`
  form is optional for checkers, and members built dynamically become attribute errors.

## Tests

- Register every marker, and keep `strict_markers` (or `strict`) on.
- A test that needs a service carries a marker. In the CI job that provides the service, a
  missing service is a failure, not a skip.
- Where the project turns warnings into errors (`filterwarnings = ["error"]`) or measures branch
  coverage, keep both on; relaxing either to make a change pass hides what the change broke.

## Environment

- After switching branches, a directory holding only `__pycache__` can import as a namespace
  package. Delete stale `__pycache__` directories, or set `PYTHONPYCACHEPREFIX` to keep caches out
  of the tree.
- Never commit `__pycache__`, virtual environments, or `.env` files.
