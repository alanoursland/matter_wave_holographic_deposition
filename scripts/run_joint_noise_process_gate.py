"""Run T52 joint dephasing, column-aberration, and process-yield gate."""

from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from t52_joint_noise_process_gate import main


if __name__ == "__main__":
    raise SystemExit(main())
