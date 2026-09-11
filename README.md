# two-ph-protein-electrostatics

Harbor tasks for [Terminal-Bench Science](https://www.terminal-bench-science.ai/): treat **pH as a first-class state variable** in protein agent workflows.

Frontier agents fold, dock, and rank proteins at a single default pH (7 or 7.4). Biology is compartmentalized. Histidine occupancy flips between endosome (pH 5.5) and cytosol (pH 7.4). These tasks make that failure numeric on **public peptides only**.

Apache-2.0. No unpublished ligands. This repository is not MoleculoSphere and is not a therapeutic program.

## Tasks

| Folder | What a pH-blind agent does | Oracle |
|---|---|---|
| `tasks/physical-sciences/chemistry/agent-debug-wrong-protonation` | ships a seed script frozen at pH 7.4 | patched code + occupancy at **both** 5.5 and 7.4 |
| `tasks/physical-sciences/chemistry/his-protonation-microstates` | one table, no second pH | `occupancy.csv` for HHHHHH and KSRRRAR |

## Bot team

Operating charter: [`BOT_TEAM.md`](BOT_TEAM.md).

## Approach to TB-Science

1. Discord `#tb-science` then `#tb-science-task-ideas`
2. [Task proposal form](https://airtable.com/appzZC5gEHrXSfNNw/pagjgS95lAQ5FVJxt/form)
3. Short email to `stevendi@stanford.edu` with this repo + the form id
4. PR only after the proposal is approved

Draft copy lives in [`authoring/`](authoring/).
