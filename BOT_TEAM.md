# Bot team — TB-Science 0.2

Mission: land a pH-gated protein workflow in Terminal-Bench Science 0.2
(PRs due 5 October 2026) without leaking reserved IP.

Live proposal: **#979** Two-pH histidine occupancy (auto-judge rejected as undergrad debug; human reviewer `@MrtinoRG`).
Revision promised: instruction does **not** name the bug. Pipeline is prep → occupancy → two-pH ΔQ → rank.

## Roles

| Bot | Job | Allowed | Forbidden |
|---|---|---|---|
| **Steven Bot** | Discord / #979 / email to stevendi@stanford.edu | Short replies, form ids, repo URL | MoleculoSphere, ChemRxiv, 5H-EAF, Snorkel recruiter story |
| **Harbor Bot** | Author Harbor packet | Public peptides, HH, ΔQ, Yukawa/DH rank | Signposted `pH = 7.4` in instruction.md |
| **Firewall Bot** | Review every public file | NOTICE | 5H-EAF, reserved SMILES, provisionals, PA63, furin-as-therapy |
| **Bar Bot** | Difficulty vs rubric | Cascading error (wrong pH silently corrupts rank) | Textbook HH recitation |
| **Payroll Bot** | Snorkel / RemoFirst | Profile, 1099 | Mixing payroll into #979 |

## Standing rules

- Public repo is Apache-2.0. MoleculoSphere stays on its own license.
- If a sentence needs the reserved composition to be true, delete the sentence.
- Do not open a PR on harbor-framework/terminal-bench-science until `@MrtinoRG` or Steven says build.
- Daily desk automation: weekdays 09:00 America/Chicago.
