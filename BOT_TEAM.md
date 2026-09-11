# Bot team — TB-Science / Steven Dillmann

Mission: get a pH-gated protein workflow into Terminal-Bench Science 0.2
(PRs due 5 October 2026) without leaking reserved IP.

## Roles

| Bot | Job | Allowed outputs | Forbidden |
|---|---|---|---|
| **Steven Bot** | Discord, Airtable, email to stevendi@stanford.edu | 8-line pitch, form answers, 12-line email | MoleculoSphere, ChemRxiv PDF, 5H-EAF, Snorkel recruiter story |
| **Harbor Bot** | Author `instruction.md`, `task.toml`, oracle, pytest | Public peptides, HH occupancy, two-pH ΔQ | Homework-only tasks with no debug/rank step |
| **Firewall Bot** | Review every file before it is public | NOTICE reminders | 5H-EAF, SMILES of reserved ligands, provisionals, PA63, furin-as-therapy |
| **Bar Bot** | Difficulty vs TB-Science rubric | Research workflow a scientist is paid to do | Textbook HH recitation as the only task |
| **Payroll Bot** | Snorkel / RemoFirst | Profile tags, 1099 carve-out | Mixing Snorkel cash story into the Steven email |

## This week's order

1. Discord join + intro + `#tb-science-task-ideas` pitch (`authoring/discord.md`)
2. Airtable proposal (`authoring/proposal.md`)
3. This public repo as the only attachment
4. Email Steven after the form has an id (`authoring/email-steven.md`)
5. Monday/Tuesday TB-Science meeting if invited

## Standing rules

- Public repo is Apache-2.0. MoleculoSphere stays on its own source-available license.
- If a sentence needs the reserved composition to be true, delete the sentence.
- Task 2 (hardcoded pH = 7.4) is the hook. Task 1 is the prerequisite.
- Next build if reviewers call Task 1 a toy: two-pH rank flip on the same public peptides.
