# two-ph-yukawa-pdb-rank

Harbor task for [Terminal-Bench Science](https://www.terminal-bench-science.ai/) discussion [#1815](https://github.com/harbor-framework/terminal-bench-science/discussions/1815).

Parse a deposited public PDB, occupy at endosome 5.5 and cytosol 7.4 from a locked table, compute intermolecular Yukawa energy versus a public acidic partner, and rank by ΔU.

Apache-2.0. Public coordinates only (`1LYZ`, `KSRRRAR`, `EEEEEE`). Not MoleculoSphere. Not a therapeutic program.

## Difficulty

The agent has to prepare a real ATOM/HETATM file (water, ligands, altLoc, missing charge centers) and then propagate two-compartment occupancy into a pair sum. A fluent paragraph about pH-dependent binding plus a single default charge still writes a plausible energy table; the numeric gate fails unless occupancy actually moves and is used in U.

## Reference solution

Stdlib PDB parse with the job-file prep rules, Henderson–Hasselbalch occupancy from the locked table, intermolecular Yukawa pair sum, rank by ΔU.

## Verification

Independent pytest recomputes occupancy, U, ΔU, and rank from `SPEC.md`. Binary reward. The ranking rule is in `SPEC.md`; the expected rank is not.
