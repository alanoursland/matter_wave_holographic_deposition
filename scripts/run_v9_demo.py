"""Run the historical v9 patterned-substrate configuration."""

from pathlib import Path
import sys

SRC = Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from iqs.pipelines.patterned_substrate import main


if __name__ == "__main__":
    main()
