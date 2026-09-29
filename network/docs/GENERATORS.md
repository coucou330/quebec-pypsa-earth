# Generators and storage

Fleet composition on the main working network (`elec_solved.nc`, 205 buses). Counts are real
units after `reduce_voltage_network.py`'s reassignment (none dropped); capacities are nameplate.

| Carrier | Count | Capacity (MW) | Source |
|---|---|---|---|
| Run-of-river hydro (`ror`) | 35 | 13,570 | Real HQ hydro station list + 2022 hourly `Hydraulique` dispatch |
| Hydro storage (`StorageUnit`, `hydro`) | 22 | 28,583 | Same station list; reservoir plants modeled with storage |
| Onshore wind (`onwind`) | 33 | 3,933 | Real HQ wind farm list (44 IPP-operated farms, official total); dispatched by weather-derived capacity factor |
| Solar (`solar`) | 3 | 10 | Real facility list (2 official stations: Gabrielle-Bodis 8 MW, Robert-A.-Boyd 2 MW); dispatched by weather-derived capacity factor |
| OCGT (gas) | 1 | 411 | Real 2022 hourly `Thermique` dispatch |
| Load shedding | 189 | (sized to local peak) | Synthetic VOLL placeholder, $10,000/MWh, one per load bus |
| Slack placeholder (`AC`) | 1 | 7,722 | Zero-dispatch generator so largest hydro unit is assigned as slack, Pypsa do not let storage unit be assigned slack bus |

Load-shedding generators are a modeling device, not real capacity: they exist so LOPF can shed
load at a heavy cost penalty instead of failing to solve, and their dispatch is the metric used to
report unserved load (currently 0% on the solved network).

Churchill Falls' real 5,428 MW is included in the hydro storage total above (`465
hydro-Churchill-Falls`, at bus `465-735kv`, added by `add_churchill_falls_tie.py`).

**Hydro capacity correction:** the original attach was missing 20 real HQ stations (~1,627 MW)
that don't name-match between `hq_major_facilities_2023.csv` and `hq_hydro_stations_official.csv`
-- either absent from the major-facilities list entirely, or present there with no coordinates.
14 of them (1,542 MW, including Eastmain-1 at 480 MW) were added manually, each attached to the
bus of an already-modeled station sharing the same river (e.g. Eastmain-1 -> Bernard-Landry's bus,
both on the Eastmain river) -- no plant-specific inflow data exists for them, so run-of-river units
are left at the default `p_max_pu=1.0` and reservoir units at zero inflow (matching
`attach_real_generators.py`'s own no-inflow fallback). The remaining 6 (~160 MW: Chute-Hemmings,
Drummondville, Lac-Robertson, Mitis-1, Mitis-2, Sept-Chutes) have no modeled station on their
river to anchor to and were left out. Hydro capacity now totals 42,153 MW against HQ's official
42,313 MW (36,885 MW HQ-owned stations + 5,428 MW Churchill Falls, excluding 706 MW of
IPP-operated hydro that isn't modeled at all).

**Wind and solar capacity correction:** both were scaled uniformly (per-generator, within each
carrier) to match HQ's official installed totals -- wind from 3,721.8 MW (39 units, pre-island-
extraction) to 3,933 MW (44 IPP-operated farms), solar from 12.3 MW (3 units) to 10.0 MW (2
stations). Station counts weren't reconciled 1:1 against the official lists, only the aggregate
capacity.

## Dispatch ceilings

Hydro (ror) and gas (OCGT) are bounded by their own real 2022 hourly dispatch from Hydro-Québec's
published generation-by-source data, not by a flat capacity factor -- `p_max_pu` at each hour is
that unit's share of the fleet-wide real dispatch ratio for its carrier at that hour, with a 10%
margin on top (`CEILING_MARGIN`) so the model isn't forced to never exceed history. Wind and solar
use standard weather-derived capacity factors from PyPSA-Earth's own renewable resource pipeline
(unchanged from the upstream default). Hydro storage reservoirs are bounded by the same real
fleet-wide ratio recovered from the hourly inflow data already attached to each unit.

## Marginal costs

Assigned per carrier from PyPSA-Earth's default `resources/costs_2030_elec.csv` (2030 cost
projections -- not Quebec-specific, used only to rank dispatch order between carriers). OCGT is
set to a small negative marginal cost so the solver prefers dispatching it up to its real ceiling
over shedding load; onshore wind is set to zero.

## Slack bus

Assigned to whichever bus holds the largest *total* generation capacity (Generator + StorageUnit
combined) on the network's largest AC-connected component. PyPSA's own bus-control logic
(`find_bus_controls()` / `find_slack_bus()`) only ever reads `Generator.control`, never
`StorageUnit.control` -- confirmed directly from `pypsa/pf.py` source -- so a storage-only bus
gets a zero-dispatch placeholder `Generator` added so it can still be a slack candidate on equal
terms. On the current solved network this lands on bus 339 (7,722 MW combined: the
`339 hydro-La-Grande-2-A` and `339 hydro-Robert-Bourassa` storage units), via the placeholder
`339 slack-placeholder`.
