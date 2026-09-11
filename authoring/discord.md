# Discord pitch (paste into #tb-science-task-ideas)

Hi — Esteban Fridman, MD-PhD (Weill Cornell neuroscience; current work is physical chemistry of pH-dependent proteins).

Proposal: science agents treat proteins as pH-invariant. They default to pH 7.4 and skip histidine occupancy between endosome (5.5) and cytosol (7.4). TB-Science 0.1 has no titration / two-pH occupancy task.

Draft Harbor packets (Apache-2.0, public peptides HHHHHH and KSRRRAR only):
https://github.com/eafridman-afk/two-ph-protein-electrostatics

Lead task is a seed script that hardcodes pH = 7.4. The oracle fails unless the agent uses the job pH.
