"""Compatibility wrapper for :mod:`iqs.sources.coherent_matterwave`.

New code should import beam models from ``iqs.sources``.  This module is
retained so historical studies and external notebooks continue to run.
"""

from iqs.sources.coherent_matterwave import *  # noqa: F401,F403


if __name__ == "__main__":
    import runpy
    from pathlib import Path

    runpy.run_path(
        str(Path(__file__).parent / "iqs" / "sources" / "coherent_matterwave.py"),
        run_name="__main__",
    )
