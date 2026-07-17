"""Compatibility wrapper for :mod:`iqs.pipelines.holographic_caging`."""

from iqs.pipelines.holographic_caging import *  # noqa: F401,F403
from iqs.pipelines.holographic_caging import main as _main


if __name__ == "__main__":
    _main()
