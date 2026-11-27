# Addendum to the multi-size ED analysis plan (26 November 2026)

This note does **not** change `PREREGISTRATION.md` (sha256
`504a9e21047076de4039610f67c502f727b49b69100a74ecd681de44d8f64564`), which stays locked
and unread-edited. It records a rule for one planned point that was left out after other
results were already visible, and the outcome of that rule.

The plan lists an 18-site (18t) run at n = 0.667 (nup = 6) as feasible under the memory
budget. That point was started twice and aborted both times by the 1.2 GB memory watchdog
(25 Nov 2026, 3:21 PM ET; 26 Nov 2026, 2:42 AM ET), on a shared machine, after the other
18t/20t fillings and the locked 8–16-site extrapolation had already been seen. A third,
gated retry was started on 26 Nov 2026, 9:06 AM ET.

**Committed rule, written before the retry finished:** if the retry completes successfully,
the 18t n = 0.667 point is included. Results will be reported **both with and without it**.
The point will not be included or excluded according to how it moves the numbers.

**Outcome.** The retry completed successfully on 26 Nov 2026, 12:48:57 PM ET
(`runs/18t_nup6_U8_k.json`: status ok, E0 = −18.223899966503254, P_org = 0.02158359261849439,
deg = 1, runtime_s = 23246.5 cumulative, peak_rss_gb = 4.022). The point is therefore
included in the locked extrapolation and post-hoc controls written on 26 Nov 2026, 1:25–1:26 PM ET.
The earlier withheld run (26 Nov 2026, 2:53 AM ET), which left the point out, is kept at
`results_without_18t_nup6_20261126-0253/` (and in the package under
`ed_multisize/results_without_18t_nup6/`).

With the point included, 18t covers n ≤ 0.667 (four sizes at n = 0.611 and n = 0.667).
18t n ≥ 0.778 and 20t n ≥ 0.6 remain skipped by the pre-registered 5 GB budget, so fits at
n ≥ 0.722 still use only the already-public 10–16-site clusters. Only n ≤ 0.556 was blind.
The n = 0.611 and n = 0.667 rows are no longer provisional; they are reported both with and
without the new point in the write-up. Including the late point lowered both positive
intercepts (0.215 → 0.095 at n = 0.611, 0.109 → 0.061 at n = 0.667). That is the direction
that makes the locked numbers smaller, not larger; the point was not included because it
helped a favourite result.

**Disclosure (added 26 Nov 2026, afternoon, after the result was known).** The rule text was
first recorded in an earlier version of this file (sha256
`2070a11e90c3103ce01d6e3a329a1aafb21d3d5b58139f2ae748a24f05b39758`, project LOG.md entry
26 Nov 2026, 9:37 AM ET), while the retry was still running. That version was overwritten
when the Outcome section was added; the file's mtime after that save was 1:28:45 PM ET,
about 40 minutes after the retry ended at 12:48:57 PM ET. This disclosure was added later
the same afternoon in draft v2.3. What was preserved: the inclusion rule itself (include
the point if it completed, and report results both with and without it). That rule was set
before the result, and it is mechanical — it follows the locked plan's list of feasible
points plus "report both ways," and leaves no room to include or drop the point by how
it moved the numbers. What changed in the 1:28 PM save: the Outcome section and the
measured values were added, and "provisional / retry in progress" language was removed.
No copy of the 9:37 AM file contents was found (no git history, no backup of the file
itself). The earlier hash is recorded in LOG.md and was listed in the v2.1 package
manifest recorded there (SHA256SUMS `d5fa9870…`, 177 files). The first-run (without-18t)
tables are preserved unaltered in `ed_multisize/results_without_18t_nup6/` (hashes match
the LOG.md entry of 2:53 AM ET). Local timestamps and LOG.md are not independently
verifiable.
