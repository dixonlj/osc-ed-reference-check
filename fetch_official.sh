#!/bin/bash
# Fetch official OSC files this package compares against or wraps (NOT bundled; third-party).
# The three data files are CC-BY-NC 4.0. osc_instrument.py is official OSC instrument code.
# Pinned to dataset revision 346f184d3eff51ebf87ffc5da20c6ec74e4cd67e and checked by sha256.
# Optional local cache (same hashes): $OSC_OFFICIAL_CACHE or ../ed_reference_official_cache, if present.
set -euo pipefail
cd "$(dirname "$0")"; mkdir -p third_party
REV=346f184d3eff51ebf87ffc5da20c6ec74e4cd67e
BASE=https://huggingface.co/datasets/FINAL-Bench/OSC-Superconductor/resolve/$REV
CACHE=${OSC_OFFICIAL_CACHE:-$(cd .. && pwd)/ed_reference_official_cache}
while read -r sha path; do
  out=third_party/$(basename "$path")
  if [ -f "$out" ]; then
    :
  elif [ -f "$CACHE/$(basename "$path")" ]; then
    cp -p "$CACHE/$(basename "$path")" "$out"
  else
    curl -sfL -o "$out" "$BASE/$path"
  fi
  echo "$sha  $out" | sha256sum -c -
done <<'LIST'
126892ef7a810834ba5e72eb9609e057ca0434c6b0516081964523e90f7ecde7 validation/multisize_ED_table_kuku.md
46060ba77b8fca1c294bdf1dae8119251b4ab1b7923acc0ba1b44322a85a3146 validation/comparison_4x4_ED_vs_VMC.md
dc7e52ad1996cf823fbba314765ccb30436508bff794402e2b951978ff39fcef scoring_meta.json
6b774e529cb8d50cb68a695f5ba6bdbdcda290031833ecc21533d6a1610b5236 instrument/osc_instrument.py
LIST
