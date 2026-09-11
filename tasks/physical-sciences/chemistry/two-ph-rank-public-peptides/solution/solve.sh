#!/usr/bin/env bash
set -euo pipefail
python3 - <<'PY'
from pathlib import Path
src = Path("/app/pipeline.py")
text = src.read_text()
text = text.replace("WORKING_PH = 7.4", "WORKING_PH = None")
text = text.replace(
    "frac = frac_protonated(WORKING_PH, pk)",
    "frac = frac_protonated(spec[\"pH\"], pk)",
)
src.write_text(text)
PY
python3 /app/pipeline.py
