"""Compatibility wrapper for :mod:`iqs.holography`.

New code should import holography models and helpers from ``iqs.holography``.
This historical module name remains import-compatible for old studies.
"""

from iqs.holography import *  # noqa: F401,F403
from iqs.holography.core import main as _main


if __name__ == "__main__":
    _main()
