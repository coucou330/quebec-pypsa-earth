# -*- coding: utf-8 -*-
"""
add_shunt_impedance.py

Adds real PyPSA `ShuntImpedance` components (not the earlier zero-P
Generator-with-fixed-q_set workaround) to mitigate persistent light-load
overvoltage on the 735kV network.

Why ShuntImpedance instead of a Generator: a ShuntImpedance's reactive
output scales with the square of its own bus voltage (q = v_mag_pu**2 *
b_pu, see pypsa/pf.py), matching a real shunt reactor bank's physics --
its absorption falls off automatically as voltage sags, so it can't push
a low-voltage bus lower the way a fixed-q_set generator could. It is a
static device (b is a single float, not a time series), so it only suits
a bus whose required reactive support has one consistent sign across the
week -- a genuinely fixed-bank candidate, not a switched one.

Candidate selection and sizing
-------------------------------
A "probe": a p_nom=0, control=PV generator was added at every bus that
exceeds 1.05pu at some point (v_mag_pu_set defaults to the bus's existing
1.0pu, holding it there), and a full AC PF was run across all 168
snapshots. Its q output at every hour is the reactive power an ideal,
continuously-variable SVC would need at that bus to hold 1.0pu exactly.

From that probe run (network/divergence_data/svc_probe_q_needed.csv):
- Buses whose required Q never changes sign (min*max > 0) and is always
  negative (absorbing) are fixed-reactor candidates -- their need is
  driven by the line's own charging (roughly load-independent, per the
  Ferranti mechanism), not by how loaded the system is.
- Buses whose required Q swings sign across the week (e.g. 345, 148,
  1291, 308, 312, 25, 150, ...) need supply at some hours and absorption
  at others. A fixed shunt cannot do both -- these are explicitly
  EXCLUDED here and would need a switched or SVC-type device instead
  (out of scope for this script).

Each qualifying bus's shunt is sized at its worst (most negative, i.e.
min) required Q from the probe run, converted to susceptance via
b = Q_MVAr / v_nom_kV**2 (b_pu = b * v_nom**2 = Q at v_mag_pu ~= 1, the
same convention pypsa.pf.calculate_dependent_values uses for shunts).

Corridor-collapse buses (1291, 3582, 148, 308, 312) are excluded by
construction -- their required Q swings sign and is often positive
(supply, not absorption): they need voltage SUPPORT under heavy transfer,
which is what series compensation (apply_series_compensation.py)
addresses, not a shunt reactor.

Usage
-----
    python network/add_shunt_impedance.py --network networks/elec_735kv.nc --output networks/elec_735kv_shunt.nc
"""
import argparse
import os

import pandas as pd
import pypsa

NETWORK_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(NETWORK_DIR)
PROBE_CSV = os.path.join(NETWORK_DIR, "divergence_data", "svc_probe_q_needed.csv")


def fixed_reactor_buses(probe_csv: str = PROBE_CSV) -> pd.Series:
    """Buses whose probe-run required Q never changes sign and is always
    negative (absorbing) -- fixed shunt reactor candidates -- mapped to
    their sizing Q (MVAr, the worst/most negative hour)."""
    s = pd.read_csv(probe_csv, index_col=0)
    s.index = s.index.astype(str)
    no_flip = (s["min_Q"] * s["max_Q"]) > 0
    candidates = s[no_flip & (s["median_Q"] < 0)]
    return candidates["min_Q"]


def apply_shunt_impedance(n: pypsa.Network, probe_csv: str = PROBE_CSV) -> None:
    sizing_q = fixed_reactor_buses(probe_csv)
    sizing_q = sizing_q[sizing_q.index.isin(n.buses.index)]
    v_nom = n.buses.loc[sizing_q.index, "v_nom"]
    b = sizing_q / v_nom**2  # Siemens; b_pu = b*v_nom**2 = Q(MVAr) at v_mag_pu ~= 1
    for bus, b_val, q_val in zip(sizing_q.index, b, sizing_q):
        n.add("ShuntImpedance", f"{bus} shunt-reactor", bus=bus, g=0.0, b=b_val)
    print(f"Added {len(sizing_q)} fixed ShuntImpedance reactors "
          f"(buses: {', '.join(sizing_q.index)}), sized to fully absorb "
          "each bus's worst-hour probe-derived reactive requirement.")
    print(sizing_q.rename("sizing_Q_MVAr (negative = absorbing)").to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", default=os.path.join(BASE_DIR, "networks", "elec_735kv.nc"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--probe-csv", default=PROBE_CSV)
    args = parser.parse_args()

    n = pypsa.Network(args.network)
    apply_shunt_impedance(n, probe_csv=args.probe_csv)
    n.calculate_dependent_values()
    n.export_to_netcdf(args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
