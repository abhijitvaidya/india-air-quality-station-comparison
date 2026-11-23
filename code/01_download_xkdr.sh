#!/bin/bash
# Download the XKDR Forum India air quality database (hourly CPCB continuous-station data, Parquet) into data/raw/xkdr.
# Source: https://airquality.xkdr.org (CC BY 4.0). Needs an API key from XKDR in the environment variable XKDR_KEY.
# The paper used the export of 30 September 2026 (data to 1 September 2025); a later export may differ slightly.
# Re-runnable: skips files already present with the right size.
set -u
cd "$(dirname "$0")/.."
: "${XKDR_KEY:?set XKDR_KEY to your XKDR API key}"
B=https://airquality.xkdr.org/v1/files
mkdir -p data/raw/xkdr
curl -sf -m 60 -H "Authorization: Bearer $XKDR_KEY" "$B?format=json" -o data/raw/xkdr/files.json || { echo "file list failed"; exit 1; }
python3 -c "import json; [print(f['key'], f['bytes']) for f in json.load(open('data/raw/xkdr/files.json'))['files']]" | while read key bytes; do
  out="data/raw/xkdr/${key#v1/}"; mkdir -p "$(dirname "$out")"
  if [ -f "$out" ] && [ "$(wc -c < "$out" | tr -d ' ')" = "$bytes" ]; then continue; fi
  curl -sf -m 600 --retry 3 -H "Authorization: Bearer $XKDR_KEY" "$B/$key" -o "$out.part" && mv "$out.part" "$out" || echo "FAILED $key"
done
echo "done: $(find data/raw/xkdr -name '*.parquet' | wc -l) parquet files"
