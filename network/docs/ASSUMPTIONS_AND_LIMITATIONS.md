# Assumptions and limitations

Every generic or unmeasured assumption used in this pipeline, with the reasoning behind the
specific value chosen..

## Generic assumptions (LOPF / general)

| Assumption | Value | Where | Why |
|---|---|---|---|
| Line security margin (`s_max_pu`) | 1.0 (relaxed from PyPSA-Earth's default 0.7) | `run_lopf_main_island.py` | The 0.7 default reserves 30% headroom for N-1 contingency; this model isn't doing contingency-constrained dispatch, so the margin only produced artificial shedding. |
| Transmission capacity margin (St Clair curve) | x3 on the derived thermal envelope | `run_lopf_main_island.py`, `st_clair.py` | St Clair is a planning-level analytical envelope computed from Hydro-Quebec's own line parameters, not an empirical thermal rating -- treated as a starting point, multiplied up uniformly rather than taken as a hard limit. |
| Transformer thermal capacity | Flat 100,000 MVA | `run_lopf_main_island.py` | Real per-transformer `s_nom` from the OSM extract varies widely and binds at several major junctions carrying 38,000-155,000 MVA of real line capacity -- flattened so transformers aren't an artificial pinch point. Reactance is scaled up proportionally to keep the real per-unit value (`x_pu`) fixed. |
| Demand scale-up | x1.1142 globally | `run_lopf_main_island.py` | Calibrated against real whole-January-2022 Hydro-Québec system demand (`historique-demande-electricite-quebec.csv`, mean 32,421 MW), not just the solved week alone. Regional demand rescaling matches real HQ totals only for buses matched to one of Quebec's 17 administrative regions; buses outside that match keep smaller synthetic values, which this factor also compensates for. |
| Dispatch ceiling headroom | x1.10 on real-data-derived ratios | `run_lopf_main_island.py` | Applied to ror, OCGT, and hydro storage ceilings (all derived from real 2022 hourly dispatch ratios) so the model isn't forced to never exceed the exact historical dispatch level. Not applied to wind/solar, which use weather-derived capacity factors, not dispatch history. |
| Link 4349 (329-3975 DC tie) capacity | x3 | `run_lopf_main_island.py` | Undersized at its source value, forcing shedding upstream despite real generation capacity existing to cover it. Insensitive to going further (x6 gave the same result as x3). |
| Line reactance at voltages other than 315/345/735/765kV | PyPSA-Earth's generic default type (German textbook, 50Hz, not Quebec-specific) | `elec_full.nc`'s `type` field, unmodified | Only the two voltage tiers the reduced networks keep are overridden with Hydro-Quebec's table (`apply_hq_line_characteristics.py`, 345 treated as 315-tier and 765 as 735-tier); nothing else is touched, so lines at other voltages keep the generic default. Doesn't affect the reduced/solved networks (they only keep >=315kV lines, and the only voltages left are 315/345/735/765) -- only the unreduced reference map. The HQ table also lists 69/120/161/230kV rows, transcribed in `hq_line_characteristics_by_voltage.csv` but not applied to any line -- see [DATA_SOURCES.md](DATA_SOURCES.md). |
| Series compensation degree | 70% | `apply_series_compensation.py` | Applied to 15 lines: the ten longest lines feeding the three buses (308, 1291, 312) diagnosed as AC PF voltage-collapse points, plus five lines (162, 156, 1938, 1049, 1448) carrying the network's highest transmission angles. Real EHV series compensation typically runs 30-70% on the longest corridors; 70% is toward the high end, chosen because it's what full-demand AC PF convergence needs -- not a measured Hydro-Quebec figure for these specific lines. Real series compensation also risks sub-synchronous resonance at high compensation degrees -- not modeled here (steady-state power flow only). |
| Circuit count on 5 sending-end corridors | 3 total circuits | `correct_sending_end_circuits.py` | 114-1291, 114-148, 340-457, 66-457 and 651-25 were modeled with only 1-2 circuits despite leaving major generating stations, where Hydro-Quebec builds 3-circuit corridors. Corrected to 3 total, r/x/b recomputed from Hydro-Quebec's per-km rates accordingly. |

## Generic assumptions (AC power flow only)

| Assumption | Value | Where | Why |
|---|---|---|---|
| Load power factor | 0.95 lagging | `run_pf.py` | Reactive demand isn't in the source data; a generic industry-typical value. |
| Generator power factor | 0.9 | `run_pf.py` | Generic reactive-capability assumption, tied to nameplate capacity rather than real-time dispatch (a synchronous machine can supply close to full reactive capability near zero real output). |
| Transformer X/R ratio | 30 | `run_pf.py` | Real transformer resistance isn't in the source data; a generic value for large power transformers. |
| PV bus eligibility threshold | Topology-dependent -- widened to "every real generator bus" on the 58-bus 735kV network; a stricter `p_nom >= 100 MW, load < 10% of local generation` threshold on the 208-bus 315kV network | `run_pf.py` (`--pv-min-capacity`, `--pv-max-load-ratio`) | PV buses hold voltage with effectively unlimited reactive power. That's stabilizing on a small, heavily meshed network but destabilizing with many PV buses on a larger network with long, weak radial corridors -- the two networks needed opposite settings. |
| Reactive compensation | 70% of each bus's own reactive demand, every snapshot | `run_pf.py` | A zero-real-power PQ generator per bus, sized to track local reactive demand at each hour rather than a fixed peak-sized shunt. Improves but does not fully close the AC PF convergence gap on its own; see [POWER_FLOW.md](POWER_FLOW.md). |

## Known data limitations

- **Parallel-circuit count isn't reliably encoded in the raw data.** Real multi-circuit corridors
  show up as separate line records rather than a per-line circuit-count field, so any code path
  that trusted that field instead of counting real records would misstate corridor capacity.
- **Local-to-backbone bus reassignment uses straight-line geographic distance**, not real grid
  connectivity. This has repeatedly misassigned demand to the wrong backbone bus, sometimes
  pooling load from remote, unconnected areas onto a nearby-looking bus -- a real data-quality
  risk anywhere this reduction step is used, and a known contributor to unrealistic local stress
  in downstream results. A more accurate (graph-based) alternative exists elsewhere in the
  pipeline but hasn't been adopted for this step.
- **AC power flow does not converge on the 315kV reduced network at current real demand (0/168).**
  The divergence is localized to specific buses but not yet root-caused, and a known setup gap in
  `run_pf.py` (LOPF dispatch isn't injected as `p_set` on that network) makes its earlier AC PF
  results unreliable. The 735kV network now converges fully (168/168) -- see
  [POWER_FLOW.md](POWER_FLOW.md) for the fixes -- but light-load overvoltage on it is still open.
