This folder is **not part of the deposited package**. It is filled at run time by
`../fetch_official.sh`, which copies or downloads official OSC files from a pinned
dataset revision (346f184d) and checks their sha256:

- multisize_ED_table_kuku.md (CC-BY-NC 4.0)
- comparison_4x4_ED_vs_VMC.md (CC-BY-NC 4.0)
- scoring_meta.json (CC-BY-NC 4.0)
- osc_instrument.py (official participant instrument; sha256 6b774e52…5236)

See ../CITATIONS.md for attribution. Do not commit or zip these fetched files with the package.
