# Network versions

Three networks filtered to different voltage resolution to fit different constraints and achieve different goals

**Note on `num_parallel`:** while some corridors are already circuit-corrected upstream
(`fix_parallel_circuits_v2.py`), five 735kV corridors leaving major generating stations still had
only 1-2 of their real 3 parallel circuits modeled. `correct_sending_end_circuits.py` fixes this
on the 735kV network below -- see [POWER_FLOW.md](POWER_FLOW.md). Not yet applied to the 315kV
network.

## 1. Unreduced (`elec_full.nc`)

The raw topology from PyPSA-Earth's OSM extraction, covering every voltage level, all of
Canada. 4,013 buses, 4,546 lines, 796 transformers, 47 DC links.

Since we are only interested in the Quebec network, the
`reduce_voltage_network.py` is created to extract the quebec grid. Can be visualized by running
(`visualize_real_network_map.py`).

![Unreduced Quebec-region network -- voltage levels, substations, loads, generators](../quebec_real_network_map.png)

## 2. 315kV reduced (`elec_reduced.nc` -> `elec_solved.nc`)

210 buses, 278 lines, 22 transformers, 6 links, 205 buses after largest-island extraction, 276
lines, 72 real generators (plus per-bus load-shedding placeholders), 22 storage units, 647 loads.
Generation capacity is corrected against HQ's official generating-stations list -- see
[GENERATORS.md](GENERATORS.md).

Built by `reduce_voltage_network.py`. It keeps every bus at 315kV or above, and folds every lower-voltage local bus onto its nearest bus -- no load, generator, or storage unit is dropped, only reassigned.
`run_lopf_main_island.py`extracts the single largest connected island and solves LOPF on it.

This is the main working network: the 2022 dispatch is solved on this network via LOPF. Demand:
~30,200 MW mean over the solved week (2022-01-01 to 2022-01-07, calibrated against historical
whole-January-2022 system demand -- reduce_voltage_network.py's nearest-bus reassignment isn't
fully deterministic run-to-run, so the exact figure can drift slightly between regenerations).
**0% load shed.**

This is also the network AC power flow would need to converge on for genuinely representative
contingency analysis. It does not: **0/168** under full nonlinear AC PF -- see
[POWER_FLOW.md](POWER_FLOW.md) for similar reactive compensation method applied in the 735kV which can be implemented in the future to support the 315kV network.

![315kV reduced network -- LOPF congestion, load, shedding, and hydro dispatch](../quebec_reduced_network_map.png)

## 3. 735kV backbone (`elec_735kv.nc`)

58 buses, 109 lines. A further reduction from the solved 315kV network down to just the 735/765kV
backbone built by
`reduce_to_735kv.py`. Every other bus is reassigned to its nearest 735kV-tier bus by shortest electrical path; loads landing on the same bus are
summed into one aggregate load, and generators/storage of the same carrier are merged the same
way. Note that this network is not rerun with LOPF so that dispatch is still accurate post generator aggregation.

Purpose-built for AC power flow tractability: small and heavily meshed, so it was expected to
converge more readily than the 315kV network, making it easier to diagnose AC PF divergence. At
current (100%) demand it fully converges (**168/168**), max line loading 79.6% -- see
[POWER_FLOW.md](POWER_FLOW.md) for the fixes (circuit-count correction, series compensation). Light-load overvoltage (29/58 buses exceed 1.05pu at some point) is still open.

Line parameters (r/x/b) on this network's 315/345kV and 735/765kV lines come from Hydro-Quebec's
own line-characteristics table (`network/hq_line_characteristics_by_voltage.csv`,
`apply_hq_line_characteristics.py`) -- see [DATA_SOURCES.md](DATA_SOURCES.md).

