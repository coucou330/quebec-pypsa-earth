# -*- coding: utf-8 -*-
"""
export_to_pandapower.py

Builds a native pandapower Network directly from a solved PyPSA network,
via pandapower's create_* functions -- not a round-trip through MATPOWER
(export_to_matpower.py + pandapower.converter.from_mpc). That round-trip
loses real line length (MATPOWER's branch format has no length field, so
pandapower's importer defaults every line to length_km=1 and folds the
total ohms into the per-km slot instead), and needs a units-aware pass to
recover it. Building the pandapower net directly keeps length_km, r/x per
km, and line charging (c_nf_per_km) all physically real from the start.

Same PV/PQ/slack classification as export_to_matpower.py (run_pf.py's
prepare_for_ac_pf()), so bus type matches PyPSA's own AC PF setup:
  - Slack bus  -> pp.create_ext_grid
  - PV buses   -> pp.create_gen (voltage-controlled, vm_pu=1.0)
  - PQ buses' generators/storage dispatch -> pp.create_sgen (fixed P, Q)
Loads -> pp.create_load. ShuntImpedance -> pp.create_shunt.

One snapshot only, same as export_to_matpower.py -- pandapower's net is a
static single-operating-point model, not a time series.

Usage
-----
    python network/export_to_pandapower.py --network networks/elec_735kv.nc --output networks/elec_735kv.p
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
import run_pf  # noqa: E402 -- same PV/PQ/slack classification as the live AC PF run

DEFAULT_NETWORK = os.path.join(BASE_DIR, "networks", "elec_735kv.nc")
BASE_MVA = 100.0
F_HZ = 60.0
GEN_PF_ASSUMED = 0.85  # generic reactive-capability assumption, not measured -- same as export_to_matpower.py


def dispatch_p(t_frame_p, t_frame_p_set, name, snap) -> float:
    """Real dispatch at `snap`, preferring `p_set` over `p` when both exist
    -- see export_to_matpower.py's copy of this function for why."""
    if name in t_frame_p_set.columns:
        return float(t_frame_p_set.at[snap, name])
    if name in t_frame_p.columns:
        return float(t_frame_p.at[snap, name])
    return 0.0


def build_pandapower_net(n: pypsa.Network, snap, base_mva: float = BASE_MVA):
    import pandapower as pp

    net = pp.create_empty_network(sn_mva=base_mva, f_hz=F_HZ)

    bus_idx = {}
    for b, row in n.buses.iterrows():
        bus_idx[b] = pp.create_bus(net, vn_kv=row.v_nom, name=b, geodata=(row.x, row.y),
                                    min_vm_pu=row.get("v_mag_pu_min", 0.9), max_vm_pu=row.get("v_mag_pu_max", 1.1))

    slack_gen = n.generators.index[n.generators.control == "Slack"]
    slack_bus = n.generators.at[slack_gen[0], "bus"] if len(slack_gen) else n.buses.index[0]
    pv_buses = set(n.generators[n.generators.control == "PV"].bus)
    pp.create_ext_grid(net, bus=bus_idx[slack_bus], vm_pu=1.0, va_degree=0.0, name=f"slack ({slack_bus})")

    omega = 2 * np.pi * F_HZ
    for l, row in n.lines.iterrows():
        length = row.length if row.length > 0 else 1.0
        c_nf_per_km = (row.b / omega) * 1e9 / length  # b [S] -> C [F] -> nF/km
        i_ka = row.s_nom / (np.sqrt(3) * row.v_nom)  # S(MVA) = sqrt(3)*V(kV)*I(kA)
        pp.create_line_from_parameters(
            net, from_bus=bus_idx[row.bus0], to_bus=bus_idx[row.bus1], length_km=length,
            r_ohm_per_km=row.r / length, x_ohm_per_km=row.x / length, c_nf_per_km=c_nf_per_km,
            max_i_ka=max(i_ka, 1e-3), name=l,
        )

    for t, row in n.transformers.iterrows():
        r_ohm = row.get("r", row.x / 30.0)  # same generic X/R=30 assumption as run_pf.py/export_to_matpower.py
        # vk_percent/vkr_percent (short-circuit voltage) are x_pu/r_pu on the
        # transformer's own s_nom base, expressed as a percentage.
        vk_percent = 100.0 * row.x / row.s_nom
        vkr_percent = 100.0 * r_ohm / row.s_nom
        pp.create_transformer_from_parameters(
            net, hv_bus=bus_idx[row.bus0], lv_bus=bus_idx[row.bus1], sn_mva=row.s_nom,
            vn_hv_kv=n.buses.at[row.bus0, "v_nom"], vn_lv_kv=n.buses.at[row.bus1, "v_nom"],
            vk_percent=max(vk_percent, 0.01), vkr_percent=max(vkr_percent, 0.001),
            pfe_kw=0.0, i0_percent=0.0, name=t,
        )

    real_gens = n.generators[(n.generators.carrier != "load_shedding") & (n.generators.carrier != "reactive_compensation")]
    for g, row in real_gens.iterrows():
        p = dispatch_p(n.generators_t.p, n.generators_t.p_set, g, snap)
        qmax = row.p_nom * np.tan(np.arccos(GEN_PF_ASSUMED)) if row.p_nom > 0 else 0.0
        if row.bus == slack_bus:
            continue  # already represented by the ext_grid
        elif row.bus in pv_buses and row.control == "PV":
            pp.create_gen(net, bus=bus_idx[row.bus], p_mw=p, vm_pu=1.0,
                           max_q_mvar=qmax, min_q_mvar=-qmax, name=g)
        else:
            # PQ generators still carry a fixed reactive supply in PyPSA's own
            # AC PF setup (prepare_for_ac_pf's q_set = p_nom*tan(acos(0.9)));
            # dropping it here would understate real reactive support at PQ buses.
            q = float(row.get("q_set", 0.0))
            pp.create_sgen(net, bus=bus_idx[row.bus], p_mw=p, q_mvar=q, name=g, controllable=False)

    for s, row in n.storage_units.iterrows():
        p = dispatch_p(n.storage_units_t.p, n.storage_units_t.p_set, s, snap)
        if row.bus == slack_bus:
            continue
        pp.create_sgen(net, bus=bus_idx[row.bus], p_mw=max(p, 0.0), q_mvar=0.0, name=s, controllable=False)

    load_p = n.loads_t.p_set.loc[snap]
    pf_angle = np.arccos(0.95)
    for ld, row in n.loads.iterrows():
        p = float(load_p.get(ld, 0.0))
        pp.create_load(net, bus=bus_idx[row.bus], p_mw=p, q_mvar=p * np.tan(pf_angle), name=ld)

    for sh, row in n.shunt_impedances.iterrows():
        q_mvar = row.b * n.buses.at[row.bus, "v_nom"]**2  # MVAr at v_mag_pu=1.0
        pp.create_shunt(net, bus=bus_idx[row.bus], q_mvar=q_mvar, p_mw=0.0, name=sh)

    return net


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--network", default=DEFAULT_NETWORK)
    parser.add_argument("--output", required=True, help="Output path (.p pickle, or .json/.xlsx via pandapower.to_*)")
    parser.add_argument("--snapshot-index", type=int, default=None)
    parser.add_argument("--pv-min-capacity", type=float, default=run_pf.PV_MIN_CAPACITY_MW)
    parser.add_argument("--pv-max-load-ratio", type=float, default=run_pf.PV_MAX_LOAD_RATIO)
    args = parser.parse_args()

    import pandapower as pp

    print(f"Loading: {args.network}")
    n = pypsa.Network(args.network)
    n.calculate_dependent_values()
    n.transformers["r"] = n.transformers["x"] / 30.0

    run_pf.prepare_for_ac_pf(n, pv_min_capacity=args.pv_min_capacity,
                              pv_max_load_ratio=args.pv_max_load_ratio,
                              reactive_compensation_ratio=0.0)

    if args.snapshot_index is None:
        loading = n.lines_t.p0.abs().div(n.lines.s_nom, axis=1).max(axis=1) if not n.lines_t.p0.empty else None
        snap = loading.idxmax() if loading is not None and len(loading) else n.snapshots[0]
        print(f"No --snapshot-index given: using {snap}")
    else:
        snap = n.snapshots[args.snapshot_index]

    net = build_pandapower_net(n, snap)
    print(f"pandapower net: {len(net.bus)} buses, {len(net.line)} lines, {len(net.trafo)} trafos, "
          f"{len(net.gen)} gen, {len(net.sgen)} sgen, {len(net.ext_grid)} ext_grid, {len(net.load)} loads, "
          f"{len(net.shunt)} shunts")

    ext = os.path.splitext(args.output)[1].lower()
    if ext == ".json":
        pp.to_json(net, args.output)
    elif ext in (".xlsx", ".xls"):
        pp.to_excel(net, args.output)
    else:
        pp.to_pickle(net, args.output)
    print(f"Saved: {args.output}")


if __name__ == "__main__":
    main()
