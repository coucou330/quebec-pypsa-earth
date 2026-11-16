"""
find_bad_hours.py  (v2)

Lists the snapshots (hours) where the AC power flow diverged, using a *_pf.nc file
saved by network/run_pf.py.

run_pf.py does not store the per-snapshot convergence flag in the file, so this uses a
heuristic: a diverged hour holds huge, negative or NaN bus voltages (outside 0.2-3.0 pu),
while converged hours stay in a plausible range (about 0.6-1.6 pu on this network).
v1 used 0.5-1.5 pu and wrongly flagged converged hours with voltages above 1.5 pu.

ALWAYS CHECK: the number of "diverged" hours printed here must equal
168 minus N, where N is the "N / 168 snapshots fully converged" line printed by run_pf.py.

Usage (from the repo root, inside `pixi shell`):
    python jana_tests/find_bad_hours.py jana_tests/c315_shunt_pf.nc
    python jana_tests/find_bad_hours.py jana_tests/c315_shunt_pf.nc \
        etienne_tests/results/diagnose/elec_solved_cc_sc.nc

The optional 2nd argument is the LOPF-solved network (before the AC power flow). It is
used to report the load shedding of the optimisation, which the *_pf.nc no longer holds.
"""
import sys

import pypsa

LO, HI = 0.2, 3.0  # voltage window (pu) used to tell converged from diverged hours

if not 2 <= len(sys.argv) <= 3:
    sys.exit(__doc__)

n = pypsa.Network(sys.argv[1])

v = n.buses_t.v_mag_pu
ok = ((v > LO) & (v < HI)).all(axis=1)  # NaN compares False, so NaN hours count as bad
bad = ok.index[~ok]

print(f"\n{int(ok.sum())} / {len(ok)} hours look converged, {len(bad)} look diverged")
print("CHECK: the diverged count must equal 168 - N from run_pf ('N / 168 ... converged').")

load = n.loads_t.p_set.sum(axis=1)
rank = load.rank(ascending=False).astype(int)  # 1 = highest-demand hour

# Voltages on the converged hours only
vo = v[ok]
if len(vo):
    outside = ((vo < 0.8) | (vo > 1.2)).any(axis=1)
    print(
        f"\nConverged hours: voltage min {vo.min().min():.3f}, max {vo.max().max():.3f} pu; "
        f"{int(outside.sum())} / {len(vo)} have at least one bus outside 0.8-1.2 pu"
    )
    if outside.any():
        print(
            f"  mean load of those hours: {load[vo.index[outside]].mean():.0f} MW "
            f"(all hours: {load.mean():.0f} MW)"
        )

# Load shedding of the optimisation, read from the LOPF-solved network if given
shed = None
if len(sys.argv) == 3:
    s = pypsa.Network(sys.argv[2])
    cols = [g for g in s.generators.index if "load_shedding" in str(g)]
    if cols and not s.generators_t.p.empty:
        shed = s.generators_t.p[[c for c in cols if c in s.generators_t.p.columns]].sum(axis=1)

if len(bad) == 0:
    sys.exit("\nNo diverged hour to list.")

print(
    f"\nMean total load (MW): all hours {load.mean():.0f}, "
    f"converged {load[ok].mean():.0f}, diverged {load[~ok].mean():.0f}"
)
print("\nDiverged hours:")
print(f"{'snapshot':<22}{'load MW':>10}{'load rank':>11}{'shed MW':>10}")
for t in bad:
    sh = f"{shed[t]:.0f}" if shed is not None else "n/a"
    print(f"{str(t):<22}{load[t]:>10.0f}{rank[t]:>11}{sh:>10}")
