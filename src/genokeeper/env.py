"""The environment variables the analyser reads: `GENOKEEPER_<name>`, else the name it had
before the project was renamed from genomelab, `GENOMELAB_<name>`, so a shell profile or an
installed LaunchAgent written before the rename keeps working."""

import os
from collections.abc import Mapping

PREFIX = "GENOKEEPER_"
OLD_PREFIX = "GENOMELAB_"


def getenv(name: str, environ: Mapping[str, str] = os.environ) -> str | None:
    """`$GENOKEEPER_<name>`, else `$GENOMELAB_<name>`, else None; an empty value counts as
    unset."""
    return environ.get(PREFIX + name) or environ.get(OLD_PREFIX + name) or None
