"""Compatibility wrapper for :mod:`iqs.pipelines.patterned_substrate`."""

from iqs.pipelines.patterned_substrate import *  # noqa: F401,F403
from iqs.pipelines.patterned_substrate import main as _main


if __name__ == "__main__":
    _main()
