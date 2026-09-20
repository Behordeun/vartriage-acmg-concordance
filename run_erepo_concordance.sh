#!/usr/bin/env bash
# Reproducible eRepo concordance run for the ACMG concordance analysis.
#
# Runs analysis/validate_erepo_api.py against a pinned vartriage version with
# strict (point-system) combining and SpliceAI scores. The first run makes the
# gnomAD API calls once and persists them (hits and confirmed misses) in the
# response cache; subsequent runs read from the warm cache and finish quickly.
#
# The classifier changes that this harness previously compared as separate
# branches ("baseline" vs "remediated") are all released in vartriage >= 0.18.3
# (SVI point system, per-population gnomAD, PM2 confirmed-absence, the gnomAD
# not-found cache). This script therefore runs a single arm against whatever
# vartriage is installed in the active environment; pin that version explicitly
# for a citable result.

set -euo pipefail

# Repo root = directory this script lives in (no hardcoded absolute paths).
CONC="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$CONC"

# Python interpreter: prefer the repo venv, fall back to python3 on PATH.
if [[ -x "$CONC/.venv/bin/python" ]]; then
  PY="$CONC/.venv/bin/python"
else
  PY="$(command -v python3 || command -v python)"
fi

STAMP="$(date +%Y%m%d_%H%M)"
mkdir -p logs results

# Keep the machine awake for the long first pass, adapted to the OS.
# macOS: caffeinate. Linux with systemd: systemd-inhibit. Otherwise: run plainly.
keep_awake() {
  if command -v caffeinate >/dev/null 2>&1; then
    caffeinate -i "$@"
  elif command -v systemd-inhibit >/dev/null 2>&1; then
    systemd-inhibit --what=idle:sleep --why="eRepo concordance run" "$@"
  else
    "$@"
  fi
}

# Report the vartriage version under test so the result is attributable.
"$PY" - <<'PYEOF'
import vartriage
print(f"  vartriage under test: {getattr(vartriage, '__version__', 'unknown')} ({vartriage.__file__})")
PYEOF

cache_state() {
  "$PY" - <<'PYEOF'
import sqlite3
from pathlib import Path
db = Path("data/erepo/.api_cache.db")
if not db.exists():
    print("  gnomad cache: (none yet)")
else:
    c = sqlite3.connect(f"file:{db}?mode=ro", uri=True)
    t = c.execute("select count(*) from api_cache where source='gnomad'").fetchone()[0]
    m = c.execute(
        "select count(*) from api_cache where source='gnomad' "
        "and response_json like '%not_found%'"
    ).fetchone()[0]
    print(f"  gnomad cache: total={t} misses={m} hits={t-m}")
    c.close()
PYEOF
}

echo "### PRE-RUN cache state ###"; cache_state

keep_awake "$PY" analysis/validate_erepo_api.py \
  --erepo     data/erepo/clingen_erepo.tsv \
  --revel     data/references/revel_grch38.tsv \
  --cache-db  data/erepo/.api_cache.db \
  --spliceai  data/erepo/spliceai_erepo.tsv \
  --combining strict \
  --output    "results/erepo_strict_${STAMP}.json" \
  2>&1 | tee "logs/erepo_strict_${STAMP}.log"

echo "### POST-RUN cache state ###"; cache_state

echo
echo "DONE. Output: results/erepo_strict_${STAMP}.json"
echo "Log:         logs/erepo_strict_${STAMP}.log"
