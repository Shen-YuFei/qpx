"""File-name helpers shared across QPX inputs."""

from __future__ import annotations

import re
from pathlib import PurePosixPath

# Bruker ``.d`` directories are usually shipped archived (``run.d.zip``,
# ``run.d.tar``, ``run.d.tar.gz``); DIA-NN reports the bare stem, so the archive
# suffix must be stripped together with ``.d``.
_RUN_FILE_SUFFIX = re.compile(
    r"(?i)\.(?:mzml(?:\.gz)?|mzxml|raw|d(?:\.(?:zip|tar(?:\.gz)?|tgz))?|wiff|mgf|dia)$"
)


def run_file_stem(value: str) -> str:
    """Return a run basename without its acquisition-file suffix."""
    name = PurePosixPath(str(value).strip().replace("\\", "/")).name
    return _RUN_FILE_SUFFIX.sub("", name)
