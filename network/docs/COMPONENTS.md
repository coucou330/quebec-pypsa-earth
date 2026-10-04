# How each component is modeled in PyPSA

Which PyPSA object type each real-world element maps to, which attributes carry its behavior, and
where those values come from. For *why* the network is split into three tiers see
[NETWORKS.md](NETWORKS.md); for real-data sources see [DATA_SOURCES.md](DATA_SOURCES.md); for
generic/unmeasured assumptions see [ASSUMPTIONS_AND_LIMITATIONS.md](ASSUMPTIONS_AND_LIMITATIONS.md).

## Bus

`Bus`: `v_nom`, `x`/`y` (geographic coordinates), `v_mag_pu_set`, `v_mag_pu_min`/`v_mag_pu_max`.

- `v_nom` is the real nominal voltage (315/345/735/765 kV on the reduced networks).
- `v_mag_pu_set` defaults to 1.0 pu -- what a slack or PV bus holds its voltage at.
- `v_mag_pu_min`/`max` are set to 0.95/1.05 pu by `run_pf.py`'s AC PF preparation, for checking
  results against -- `n.pf()` does not enforce them.
- Bus *type* (slack/PV/PQ) is not a `Bus` attribute -- see Generator below.

## Line

`Line`: `r`, `x` [ohm], `b` [S], `length` [km], `num_parallel`, `type`, `s_nom` [MVA], `s_max_pu`.

- `r`, `x`, `b` on every 315/345kV and 735/765kV line come from Hydro-Quebec's own
  line-characteristics table (`apply_hq_line_characteristics.py`), per circuit:
  `r, x = r_per_km, x_per_km * length / num_parallel`; `b = b_per_km * 1e-6 * length * num_parallel`.
  `type` is cleared afterward so `n.calculate_dependent_values()` doesn't overwrite them from
  PyPSA-Earth's generic line-type library.
- `num_parallel`: five 735kV corridors leaving major generating stations were undercounted (only
  1-2 of their real 3 circuits); `correct_sending_end_circuits.py` corrects this and recomputes
  r/x/b accordingly.
- **Series compensation is not a separate component.** It's `x *= (1 - 0.70)` applied directly to
  15 lines' own `x` (`apply_series_compensation.py`) -- a capacitor bank modeled as reduced line
  reactance, matching how it works physically.
- `s_nom`: from the St Clair loadability curve (`st_clair.py`) with a 3x margin, using each line's
  real length, voltage and (corrected) circuit count. Only enforced by LOPF, not by AC PF -- see
  [POWER_FLOW.md](POWER_FLOW.md).

## Transformer

`Transformer`: `r`, `x` [pu on own `s_nom`], `s_nom` [MVA].

315kV network only (the 735kV backbone has none). `s_nom` is a flat 100,000 MVA (real per-transformer
values varied too widely and bound at major junctions -- see
[ASSUMPTIONS_AND_LIMITATIONS.md](ASSUMPTIONS_AND_LIMITATIONS.md)); `r` isn't in the source data, so
`run_pf.py` sets `r = x / 30` (generic X/R ratio) before any AC PF run.

## Link

`Link`: `p_set`/`p0` [MW], `p_min_pu`/`p_max_pu`, `efficiency`.

The real HVDC and back-to-back interconnections (Churchill Falls tie aside, which is a Line/bus,
not a Link). A `Link` is **not admittance-based** -- it never enters the AC or DC power-flow
Jacobian/B-matrix. Both `n.pf()` and `n.lpf()` just read `p_set` and inject it as a fixed power at
`bus0` (extraction) and `bus1` (injection x `efficiency`), exactly like a generator/load pair; no
impedance, angle drop, or thermal limit applies to a Link in power flow. In LOPF, `p0` is instead a
decision variable bounded by `p_min_pu`/`p_max_pu * p_nom` (here -1 to 1, i.e. bidirectional).

`reduce_to_735kv.py` carries no links across into the 735kV backbone -- the 735kV network has
zero `Link` components, so none of its AC PF results include interconnection flow.

## Generator

`Generator`: `p_set(t)` [MW], `p_nom`, `q_set`, `control` (`Slack`/`PV`/`PQ`), `carrier`.

- `p_set(t)` is what `n.pf()` reads (not the LOPF-solved `p`) -- `run_pf.py`'s AC PF preparation
  copies LOPF dispatch into it.
- `control` determines bus type in the power-flow solve: one `Slack` generator (the largest-
  capacity unit; picked in `run_lopf_main_island.py`, `114 ror` for 735kV AC PF work), a subset
  marked `PV` (largest generator at each bus meeting a size/load-ratio threshold -- see
  `run_pf.py`'s `--pv-min-capacity`/`--pv-max-load-ratio`), the rest left `PQ`.
- **PQ generators still carry a fixed reactive supply**: `q_set = p_nom * tan(acos(0.9))`, a
  generic reactive-capability assumption tied to nameplate capacity, not real-time dispatch.
- **PV generators hold voltage with effectively unlimited reactive power** -- PyPSA has no
  `q_max`/`q_min` attribute on `Generator` and no `enforce_q_lims`-style option in `n.pf()`. Some
  sending buses (114, 457) supply several GVAr at exactly 1.00pu, more than a real capability
  curve would allow. This is a real, unresolved modeling gap -- see
  [POWER_FLOW.md](POWER_FLOW.md).
- Storage-only buses get a zero-dispatch placeholder `Generator` (`carrier="hydro"`, `p_nom` =
  that bus's storage capacity) added purely so the bus is PV-eligible -- PyPSA's
  `find_bus_controls()` only ever reads `Generator.control`, never `StorageUnit.control`.
- The `load_shedding` carrier is a modeling device, not a real generator -- see
  [GENERATORS.md](GENERATORS.md).

## StorageUnit

`StorageUnit`: `p_set(t)` [MW], `p_nom`, `max_hours`.

Real hydro reservoirs. `p_set(t)` is the LOPF dispatch, copied in the same way as `Generator.p_set`.
Carries no `q_set` of its own -- any reactive support at a storage-only bus comes from that bus's
placeholder `Generator` (above), not the `StorageUnit` itself.

## Load

`Load`: `p_set(t)` [MW], `q_set(t)` [MVAr].

`p_set(t)` is real demand (from Hydro-Quebec municipal/regional consumption data, calibrated to
real system-wide demand -- see [DATA_SOURCES.md](DATA_SOURCES.md)). `q_set(t)` isn't in the source
data; `run_pf.py` sets it to `p_set * tan(acos(0.95))`, a generic lagging power factor.

## ShuntImpedance

`ShuntImpedance`: `g` [S] (unused, always 0 here), `b` [S].

The PyPSA component for a physical shunt reactor/capacitor bank. Not present in the committed
networks; `add_shunt_impedance.py` can add one per bus. Its power scales with the square of its own bus voltage
(`q = v_mag_pu**2 * b_pu`, set in `pypsa/pf.py`), so it's a fixed admittance, not a controllable
injection: unlike a fixed-`q_set` generator, its absorption falls off automatically as voltage
sags. `b` has no time dimension (a single float, not a time series), so it can't switch on/off or
vary with demand the way a real switched reactor bank does.

Tried (`add_shunt_impedance.py`) on the 14 persistently-overvoltage 735kV buses whose required
reactive support never changes sign, sized from a probe run (a temporary zero-`p_nom`, `v_mag_pu`-
held-at-1.0 `PV` generator at each candidate bus, reading back its required `q` over all 168
hours). **Not adopted into the committed network** -- see [POWER_FLOW.md](POWER_FLOW.md) for the
over-absorption trade-off this ran into.
