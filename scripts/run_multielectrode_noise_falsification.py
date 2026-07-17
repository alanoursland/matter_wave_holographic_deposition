"""Run the T44 full-array thermal phase-noise falsification analysis."""

from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from t44_multielectrode_noise_falsification import main


if __name__ == "__main__":
    raise SystemExit(main())
