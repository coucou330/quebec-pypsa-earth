# Meeting with Prof. Wang — November 17, 2026

Progress report on the Québec 315 kV AC model (ECSE 458 capstone), in English.

- **Slides:** https://claude.ai/artifact/3cMKQN3iFZ2HPqv4wSdnM4 (24 slides, speaker notes on each; private until shared — export to PowerPoint or PDF from the deck's menu)
- **Notebook:** `Wang_meeting_2026-11-17.ipynb` — the evidence behind every number in the slides (already executed; re-run with `jupyter nbconvert --to notebook --execute`)

## Files

| File | Content |
|---|---|
| `Wang_meeting_2026-11-17.ipynb` | the notebook (10 sections + appendix) |
| `make_results.py` | builds v2 (v1 + run-of-river priority) and writes everything in `data/` (~3 min; needs `../reseau_315kv/reseau315.py` and Nathan's `network/networks_current/*.nc`, read only) |
| `figs.py` | all figures of the notebook and of the slides (`figures/`) |
| `_make_notebook.py` | generates the notebook from text (edit here, then re-run) |
| `data/quebec_315kv_v2_peak.m` | MATPOWER case of the v2 peak hour (stable), for AC-CFM |
| `data/reseau_315kv_v2.nc` | the v2 network (PyPSA) |
| `data/*.csv` | week results, voltages, regulators, P-V curves, validation, sensitivities, below-315 kV tests |

## Definitions

- **v1** = corrected model of the Nov. 16-17 night (`../reseau_315kv/`).
- **v2** = v1 + run-of-river plants at −0.2 $/MWh in the LOPF (recommended base case).
- **Converged** = Newton converges AND every machine/SVC within its Q limits AND every bus within 0.90-1.10 pu AND voltage-stable (dV/dQ > 0 at every bus).

## AI use

Most code and diagnostics were produced with Claude (Anthropic) under Jean's direction and reviewed by a second model acting as a critical reviewer (`../revue_fable/RAPPORT.md`). The team has not yet re-derived every step by hand: that is the next phase (notebook, Section 10). Nothing was written in `network/networks_current/`.
