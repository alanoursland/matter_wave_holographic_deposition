"""Run the T49 correlated process-yield gate under quantum dephasing."""

from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from t49_dephased_process_yield_gate import main


if __name__ == "__main__":
    raise SystemExit(main())
