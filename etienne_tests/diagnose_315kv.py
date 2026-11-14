#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
diagnose_315kv.py

One-shot diagnostic for the 315 kV network's AC power-flow convergence.

STAGE 1 "baseline"  -- elec_solved.nc as it is today:
    * DC power flow (lpf): worst-angle lines, count of lines > 30 deg
    * AC power flow (pf) under two PV-bus settings (run_pf defaults, and
      all-generators-PV = the flags used for the 735 kV network):
        - how many of the 168 snapshots converge
        - voltage range over CONVERGED snapshots only (diverged snapshots hold
          garbage like 1e29 pu and would poison any summary)
        - which snapshots / hours / sub-networks fail

STAGE 2 "ported"    -- the two 735 kV fixes applied to the 315 kV network:
    elec_reduced.nc -> correct_sending_end_circuits.py
                    -> apply_series_compensation.py
                    -> run_lopf_main_island.py          (re-optimise: x changed)
    and then the same reports as stage 1.

At the end it prints one comparison table (stage x PV setting).

Run from anywhere (it locates the repo root itself):

    python etienne_tests/diagnose_315kv.py
    python etienne_tests/diagnose_315kv.py --baseline-only
    python etienne_tests/diagnose_315kv.py --ported-only --solver highs

Outputs go to etienne_tests/results/diagnose/ (CSV tables + intermediate .nc).
Nothing in network/networks_current/ is modified.

Needs PyPSA 0.30.x and the repo's own network/*.py scripts.
"""
import argparse
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pypsa


def find_repo_root():
    here = os.path.dirname(os.path.abspath(__file__))
    d = here
    for _ in range(5):
        if os.path.isdir(os.path.join(d, "network")) and os.path.isfile(
            os.path.join(d, "network", "run_pf.py")
        ):
            return d
        d = os.path.dirname(d)
    raise SystemExit("Could not find the repo root (a folder containing network/run_pf.py).")


REPO = find_repo_root()
NETWORK_DIR = os.path.join(REPO, "network")
sys.path.insert(0, NETWORK_DIR)

import run_pf  # noqa: E402  (the repo's own AC-PF preparation code)

ANGLE_FLAG_DEG = 30.0
SANE_V_RANGE = (0.8, 1.2)  # a "converged" snapshot outside this is treated as suspect

PV_VARIANTS = [
    (run_pf.PV_MIN_CAPACITY_MW, run_pf.PV_MAX_LOAD_RATIO, "run_pf defaults"),
    (0.0, 1.0, "all gens PV (735kV flags)"),
]


# --------------------------------------------------------------------------
# helpers
# --------------------------------------------------------------------------
def load_fix_lists():
    """Line IDs already covered by the 735 kV circuit-count / compensation fixes."""
    circ, comp = set(), set()
    try:
        import correct_sending_end_circuits as c

        circ = {str(k) for k in getattr(c, "NEW_NUM_PARALLEL", {})}
    except Exception as e:  # noqa: BLE001
        print(f"[note] could not read NEW_NUM_PARALLEL: {e}")
    try:
        import apply_series_compensation as a

        comp = {str(k) for k in getattr(a, "TARGET_LINES", {})}
    except Exception as e:  # noqa: BLE001
        print(f"[note] could not read TARGET_LINES: {e}")
    return circ, comp


def add_missing_slacks(m):
    """PyPSA crashes with KeyError: None if a sub-network has no generator at all.
    Give each such sub-network an inert p_nom=0 slack placeholder."""
    m.determine_network_topology()
    added = []
    for sn in m.sub_networks.index:
        buses = m.buses.index[m.buses.sub_network == sn]
        if not m.generators.bus.isin(buses).any():
            b = buses[0]
            name = f"{b} diag-slack-placeholder"
            if name not in m.generators.index:
                m.add("Generator", name, bus=b, carrier="AC", p_nom=0.0, control="Slack")
                added.append(b)
    return added


def prepare(n, pv_min=None, pv_ratio=None, ac=True):
    """Copy of n with dispatch in p_set, ready for lpf (ac=False) or pf (ac=True)."""
    m = n.copy()
    run_pf.use_lopf_dispatch_as_setpoints(m)
    m.calculate_dependent_values()
    m.determine_network_topology()
    if ac:
        run_pf.prepare_for_ac_pf(m, pv_min_capacity=pv_min, pv_max_load_ratio=pv_ratio)
    added = add_missing_slacks(m)
    if added:
        print(f"  [fix] slack placeholders added for generator-less sub-networks at buses: {added}")
    return m


def sh(cmd):
    print("\n$ " + " ".join(cmd))
    subprocess.run(cmd, check=True, cwd=REPO)


# --------------------------------------------------------------------------
# DC report: angles and loading
# --------------------------------------------------------------------------
def dc_report(n, label, top, outdir, fixes):
    m = prepare(n, ac=False)
    m.lpf(m.snapshots)

    va = m.buses_t.v_ang
    d = np.degrees(va[m.lines.bus0].values - va[m.lines.bus1].values)
    ang = pd.DataFrame(np.abs(d), index=m.snapshots, columns=m.lines.index)
    loading = m.lines_t.p0.abs().div(m.lines.s_nom, axis=1)

    circ, comp = fixes
    T = pd.DataFrame(
        {
            "bus0": m.lines.bus0,
            "bus1": m.lines.bus1,
            "v_nom": m.lines.v_nom,
            "length_km": m.lines.length.round(1),
            "circuits": m.lines.num_parallel,
            "max_angle_deg": ang.max().round(1),
            "max_loading": loading.max().round(3),
        }
    )
    T["km_per_circuit"] = (T.length_km / T.circuits).round(1)
    T["in_735kV_fix_lists"] = [
        ",".join(tag for tag, s in (("circuits", circ), ("compensation", comp)) if str(i) in s)
        for i in T.index
    ]
    T = T.sort_values("max_angle_deg", ascending=False)

    n_over = int((T.max_angle_deg > ANGLE_FLAG_DEG).sum())
    print(f"\n--- DC power flow [{label}] ---")
    print(f"  lines over {ANGLE_FLAG_DEG:.0f} deg: {n_over} / {len(T)};  "
          f"largest angle {T.max_angle_deg.iloc[0]:.1f} deg (line {T.index[0]});  "
          f"max loading {T.max_loading.max():.1%}")
    print(f"\n  Top {top} lines by max angle "
          f"(blank last column = NOT covered by the 735 kV fix lists -> candidates to inspect):")
    print(T.head(top).to_string())

    safe = label.replace(" ", "_")
    T.to_csv(os.path.join(outdir, f"worst_angles_{safe}.csv"))
    return n_over, float(T.max_angle_deg.iloc[0])


# --------------------------------------------------------------------------
# AC report: convergence
# --------------------------------------------------------------------------
def ac_report(n, label, pv_min, pv_ratio, pv_tag, outdir):
    print(f"\n--- AC power flow [{label}] PV setting: {pv_tag} "
          f"(min cap {pv_min:.0f} MW, load ratio {pv_ratio:.0%}) ---")
    m = prepare(n, pv_min, pv_ratio, ac=True)
    m.lpf(m.snapshots)  # DC seed, as run_pf.py does
    res = m.pf(m.snapshots, use_seed=True)

    conv = res.converged.all(axis=1)
    n_conv, n_all = int(conv.sum()), len(conv)

    v = m.buses_t.v_mag_pu.loc[conv]
    sane = pd.Series(False, index=conv.index)
    if len(v):
        sane.loc[v.index] = ((v > SANE_V_RANGE[0]) & (v < SANE_V_RANGE[1])).all(axis=1)
    n_sane = int(sane.sum())

    print(f"  converged: {n_conv} / {n_all}   "
          f"(of those, {n_sane} have every bus within {SANE_V_RANGE[0]}-{SANE_V_RANGE[1]} pu)")
    if len(v):
        vs = m.buses_t.v_mag_pu.loc[sane]
        vsel = vs if len(vs) else v
        print(f"  voltage over {'sane' if len(vs) else 'converged'} snapshots: "
              f"min {vsel.min().min():.3f}, max {vsel.max().max():.3f} pu;  "
              f"buses ever <0.95: {int((vsel.min() < 0.95).sum())}, "
              f">1.05: {int((vsel.max() > 1.05).sum())} of {vsel.shape[1]}")
        vmin, vmax = float(vsel.min().min()), float(vsel.max().max())
    else:
        print("  no converged snapshots -> no voltage statistics")
        vmin = vmax = float("nan")

    bad = conv.index[~conv]
    if len(bad):
        print(f"  failing snapshots: {len(bad)}  "
              f"(first {bad.min()}, last {bad.max()})")
        hours = getattr(bad, "hour", None)
        if hours is not None:
            hc = pd.Series(hours).value_counts().sort_index()
            print("  failures by hour of day: " + ", ".join(f"{h:02d}h:{c}" for h, c in hc.items()))
        tot = m.loads_t.p_set.sum(axis=1)
        print(f"  mean total load (MW): failing {tot.loc[bad].mean():,.0f}  vs  "
              f"converged {tot.loc[conv].mean() if n_conv else float('nan'):,.0f}")
        fails = (~res.converged).sum()
        for sn, k in fails[fails > 0].items():
            nb = int((m.buses.sub_network == sn).sum())
            print(f"  sub-network {sn}: {k} failed snapshots ({nb} buses)")
    safe = f"{label}_{pv_tag}".replace(" ", "_").replace("(", "").replace(")", "")
    pd.DataFrame({"converged": conv, "sane": sane}).to_csv(
        os.path.join(outdir, f"convergence_{safe}.csv")
    )
    return {
        "converged": n_conv,
        "sane": n_sane,
        "v_min": round(vmin, 3),
        "v_max": round(vmax, 3),
    }


def report_stage(label, path, top, outdir, fixes, summary):
    print("\n" + "=" * 78)
    print(f"STAGE: {label}\n  network: {path}")
    print("=" * 78)
    n = pypsa.Network(path)
    print(f"  {len(n.buses)} buses, {len(n.lines)} lines, {len(n.snapshots)} snapshots")
    n_over, max_ang = dc_report(n, label, top, outdir, fixes)
    for pv_min, pv_ratio, tag in PV_VARIANTS:
        r = ac_report(n, label, pv_min, pv_ratio, tag, outdir)
        summary.append(
            {"stage": label, "PV setting": tag, **r, "lines>30deg": n_over,
             "max_angle_deg": round(max_ang, 1)}
        )


# --------------------------------------------------------------------------
# stage 2: port the 735 kV fixes to the 315 kV network
# --------------------------------------------------------------------------
def build_ported(args):
    os.makedirs(args.outdir, exist_ok=True)
    cc = os.path.join(args.outdir, "elec_reduced_cc.nc")
    sc = os.path.join(args.outdir, "elec_reduced_cc_sc.nc")
    solved = os.path.join(args.outdir, "elec_solved_cc_sc.nc")
    py = sys.executable
    print("\nBuilding the 'ported' network (circuit fix -> series compensation -> LOPF).")
    print("Watch for '[warn] target line(s) not present': those hardcoded 735 kV line IDs "
          "do not exist in the 315 kV network.")
    sh([py, os.path.join(NETWORK_DIR, "correct_sending_end_circuits.py"),
        "--network", args.reduced, "--output", cc])
    sh([py, os.path.join(NETWORK_DIR, "apply_series_compensation.py"),
        "--network", cc, "--output", sc])
    lopf = [py, os.path.join(NETWORK_DIR, "run_lopf_main_island.py"),
            "--network", sc, "--output", solved]
    if args.solver:
        lopf += ["--solver", args.solver]
    sh(lopf)
    return solved


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawTextHelpFormatter)
    p.add_argument("--solved", default=os.path.join(NETWORK_DIR, "networks_current", "elec_solved.nc"),
                   help="LOPF-solved 315 kV network (baseline)")
    p.add_argument("--reduced", default=os.path.join(NETWORK_DIR, "networks_current", "elec_reduced.nc"),
                   help="pre-LOPF 315 kV network (input to the ported-fixes stage)")
    p.add_argument("--outdir", default=os.path.join(REPO, "etienne_tests", "results", "diagnose"))
    p.add_argument("--top", type=int, default=20, help="how many worst-angle lines to print")
    p.add_argument("--solver", default=None, help="solver for run_lopf_main_island.py (default: its own)")
    p.add_argument("--baseline-only", action="store_true")
    p.add_argument("--ported-only", action="store_true")
    args = p.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    fixes = load_fix_lists()
    summary = []

    if not args.ported_only:
        report_stage("baseline", args.solved, args.top, args.outdir, fixes, summary)

    if not args.baseline_only:
        try:
            ported = build_ported(args)
        except subprocess.CalledProcessError as e:
            print(f"\n[error] a pipeline step failed ({e}). Stage 2 skipped; see its output above.")
            ported = None
        if ported:
            report_stage("ported fixes", ported, args.top, args.outdir, fixes, summary)

    if summary:
        S = pd.DataFrame(summary)
        print("\n" + "=" * 78)
        print("SUMMARY (convergence counted over all snapshots; voltages over sane converged ones)")
        print("=" * 78)
        print(S.to_string(index=False))
        S.to_csv(os.path.join(args.outdir, "summary.csv"), index=False)
        print(f"\nCSV tables written to {args.outdir}")


if __name__ == "__main__":
    main()