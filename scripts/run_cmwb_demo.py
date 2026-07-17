"""Run the coherent matter-wave beam demonstrations."""

import runpy
from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))


if __name__ == "__main__":
    runpy.run_module("iqs.sources.coherent_matterwave", run_name="__main__")
