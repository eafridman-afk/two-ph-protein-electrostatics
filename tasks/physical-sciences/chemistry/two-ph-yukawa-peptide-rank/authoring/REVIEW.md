# Review — SPEC.md + instruction.md (1LYZ Harbor packet)

Reviewed 16 September 2026 against the locked 1LYZ approach (not 2LZM, not α7, not furin, not ATP7A).

## Keep

- Ranking **rule** only. Neither file names a winner.
- Raw deposited `1LYZ.pdb`, ignore HETATM/HOH, chain A, altloc blank/`A`.
- Charge centers HIS-ND1, LYS-NZ, ARG-CZ, GLU-CD.
- Missing named atom → q = 0, `missing_center=true`.
- Locked table occupancy; ionic strength enters κ only.
- Intermolecular Yukawa vs public `EEEEEE`. Candidates `1LYZ` and `KSRRRAR`.
- Zhao 2023 cited as a design rule from the abstract, not as a delivery claim.
- Four artifacts, binary pytest, no PROPKA/APBS/web.

## Fixes applied before freeze

1. **ASP removed from the contract.** Draft SPEC listed ASP-CG and pKa 4.00. The locked table is His 6.00, Lys 10.40, Arg 12.50, Glu 4.40 only. Asp occupancy would move in the 5.5/7.4 window and muddy the His15 term. Asp is skipped (no row).
2. **Occupancy rows = titratable residues only** (HIS, LYS, ARG, GLU), as the draft said “all others skip.” Oracle and tests match that, not 129 protein rows.
3. **`job.json` carries `candidates` and `chains`** so instruction’s “specified chain” is on disk (`1LYZ` → A).
4. **No rank-1, no HHHHHH, no 2LZM, no ATP7A, no receptor.**

## Residual (acceptable)

- `energy.csv` / `rank.csv` still use column name `peptide` for `1LYZ`. Ugly, matches proposal #994.
- Zhao sentence states the occupancy idea. It does not name the numeric winner.

Oracle-only (not in SPEC or instruction): His15 occupancy is the term that moves ΔU vs `EEEEEE` between 5.5 and 7.4 relative to the Arg-locked control.
