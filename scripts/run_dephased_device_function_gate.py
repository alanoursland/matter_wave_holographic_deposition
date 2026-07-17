"""Run the T48 electrical device gate under quantum dephasing."""

from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from t48_dephased_device_function_gate import main


if __name__ == "__main__":
    raise SystemExit(main())
