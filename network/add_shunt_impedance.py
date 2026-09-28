# -*- coding: utf-8 -*-
"""
add_shunt_impedance.py

Adds real PyPSA `ShuntImpedance` components to mitigate persistent
light-load/always-overvoltage buses on the 735kV network. A ShuntImpedance's
reactive output scales with the square of its own bus voltage (q =
v_mag_pu**2 * b_pu, see pypsa/pf.py), matching a real shunt reactor bank's
physics -- its absorption falls off automatically as voltage sags, unlike a
fixed-q_set generator. It is a static device (b is a single float, not a
time series), so it only suits a bus whose required reactive support has
one consistent sign across the week.

Candidate selection and sizing
-------------------------------
A "probe": a p_nom=0, control=PV generator is added at every bus that
exceeds 1.05pu at some point, and a full AC PF is run across all 168
snapshots. Its q output at every hour is the reactive power an ideal,
continuously-variable SVC would need at that bus to hold 1.0pu exactly.

Buses whose required Q never changes sign (min*max > 0) and is always
negative (absorbing) are fixed-reactor candidates. Buses whose required Q
swings sign across the week need supply at some hours and absorption at
others -- a fixed shunt cannot do both, so they are excluded here and would
need a switched or SVC-type device instead (out of scope for this script).

Each qualifying bus's shunt is sized at its worst (most negative) required
Q from the probe run, converted to susceptance via b = Q_MVAr / v_nom_kV**2
(b_pu = b * v_nom**2 = Q at v_mag_pu ~= 1, the same convention
pypsa.pf.calculate_dependent_values uses for shunts).

Usage
-----
    python network/add_shunt_impedance.py --network networks/elec_735kv.nc --output networks/elec_735kv_shunt.nc
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
import pypsa

NETWORK_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(NETWORK_DIR)
sys.path.insert(0, NETWORK_DIR)


def probe_reactive_requirement(n: pypsa.Network, pv_min_capacity: float = 0.0,
                                pv_max_load_ratio: float = 1.0) -> pd.DataFrame:
    """Runs a full AC PF with an ideal zero-p, 1.0pu-held SVC probe at every
    bus that exceeds 1.05pu at some point, and returns each probed bus's
    required reactive power (MVAr) over the week."""
    import run_pf

    n = n.copy()
    n.calculate_dependent_values()
    n.determine_network_topology()
    run_pf.prepare_for_ac_pf(n, pv_min_capacity=pv_min_capacity, pv_max_load_ratio=pv_max_load_ratio)
    n.lpf(n.snapshots)
    n.pf(n.snapshots, use_seed=True)

    overvoltage_buses = n.buses_t.v_mag_pu.columns[(n.buses_t.v_mag_pu > 1.05).any()]
    probes = []
    for b in overvoltage_buses:
        name = f"{b} svc-probe"
        n.add("Generator", name, bus=b, carrier="svc_probe", p_nom=0.0, control="PV")
        probes.append(name)

    n.lpf(n.snapshots)
    n.pf(n.snapshots, use_seed=True)

    q = n.generators_t.q.loc[:, probes]
    q.columns = [c.split()[0] for c in q.columns]
    return pd.DataFrame({"min_Q": q.min(), "median_Q": q.median(), "max_Q": q.max()})


def fixed_reactor_sizing(probe: pd.DataFrame) -> pd.Series:
    """Buses whose probe-run required Q never changes sign and is always
    negative (absorbing) -- fixed shunt reactor candidates -- mapped to
    their sizing Q (MVAr, the worst/most negative hour)."""
    no_flip = (probe["min_Q"] * probe["max_Q"]) > 0
    candidates = probe[no_flip & (probe["median_Q"] < 0)]
    return candidates["min_Q"]


def apply_shunt_impedance(n: pypsa.Network, sizing_q: pd.Series) -> None:
    sizing_q = sizing_q[sizing_q.index.isin(n.buses.index)]
    v_nom = n.buses.loc[sizing_q.index, "v_nom"]
    b = sizing_q / v_nom**2  # Siemens; b_pu = b*v_nom**2 = Q(MVAr) at v_mag_pu ~= 1
    for bus, b_val in zip(sizing_q.index, b):
        n.add("ShuntImpedance", f"{bus} shunt-reactor", bus=bus, g=0.0, b=b_val)
    print(f"Added {len(sizing_q)} fixed ShuntImpedance reactors "
          f"(buses: {', '.join(sizing_q.index)}), sized to fully absorb "
          "each bus's worst-hour probe-derived reactive requirement.")
    print(sizing_q.rename("sizing_Q_MVAr (negative = absorbing)").to_string())


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", default=os.path.join(BASE_DIR, "networks", "elec_735kv.nc"))
    parser.add_argument("--output", required=True)
    parser.add_argument("--pv-min-capacity", type=float, default=0.0)
    parser.add_argument("--pv-max-load-ratio", type=float, default=1.0)
    args = parser.parse_args()

    n = pypsa.Network(args.network)
    probe = probe_reactive_requirement(n, pv_min_capacity=args.pv_min_capacity,
                                        pv_max_load_ratio=args.pv_max_load_ratio)
    sizing_q = fixed_reactor_sizing(probe)
    apply_shunt_impedance(n, sizing_q)
    n.calculate_dependent_values()
    n.export_to_netcdf(args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
