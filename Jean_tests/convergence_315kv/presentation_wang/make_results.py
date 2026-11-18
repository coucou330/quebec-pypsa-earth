"""Computes everything the Nov. 17 presentation shows that is not already in ../reseau_315kv/resultats/.

v1 = "final case" of Nov. 16-17 night (../reseau_315kv/resultats/), unchanged.
v2 = v1 + ONE change found by the review (../revue_fable/RAPPORT.md, B4): run-of-river plants get a
     marginal cost of -0.2 $/MWh in the LOPF (they cannot store water), so the degenerate LOPF no longer
     shuts down Beauharnois / Carillon at the peak.
Writes data/*.csv, data/reseau_315kv_v2.nc and data/quebec_315kv_v2_peak.m. Nothing in ../reseau_315kv is modified.
Runtime ~8 min."""
import os, sys, time, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "reseau_315kv"))
import numpy as np, pandas as pd
import reseau315 as R

OUT = os.path.join(HERE, "data"); os.makedirs(OUT, exist_ok=True)
PEAK = pd.Timestamp("2022-01-03 17:00")
t0 = time.time()

def log(*a):
    print(f"[{time.time() - t0:5.0f}s]", *a, flush=True)

# ------------------------------------------------------------------ v2 build
_lopf0 = R.lopf
def lopf_ror_priority(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _lopf0(n, solver=solver, log=log)

full = R.load(os.path.join(R.CURRENT, "elec_full.nc"))
T = R.mercator_table(full)
T.to_csv(os.path.join(OUT, "mercator_table.csv"))
R.lopf = lopf_ror_priority
with R.quiet():
    n, info = R.build(T=T)
R.lopf = _lopf0
eng = R.ACEngine(n, info["ac"]["regulators"])
log("v2 built:", info["lopf"]["status"], "shed %.4f %%" % info["lopf"]["shed_pct"])

# ------------------------------------------------------------------ week
Tt, V, Q, K, C, S = R.solve_week_reserve(eng)
log("v2 week solved")

def jac(eng, k, c, Vc, st):
    Y = eng.Y0 + np.diag(1j * (k * eng.b_react + c * eng.b_cap))
    pv = sorted({eng.reg_bus[j] for j in range(len(eng.regs)) if st[j] == 0}); pvs = set(pv)
    pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]; pvpq = np.r_[pv, pq].astype(int)
    dth, dvm = eng._dSdV(Y, Vc)
    J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real],
                  [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
    return J, len(pvpq), pq

def vq_stability(eng, k, c, Vc, st):
    """V-Q sensitivity criterion: dV_i/dQ_i > 0 at every PQ bus and all eigenvalues of the reduced
    Jacobian J_R = dQ/dV positive (Gao, Morison & Kundur 1992)."""
    J, m, pq = jac(eng, k, c, Vc, st)
    Ji = np.linalg.inv(J); d = np.diag(Ji[m:, m:])
    ev = np.linalg.eigvals(np.linalg.inv(Ji[m:, m:])).real
    return int((d < 0).sum()), float(ev.min())

rows = []
for sn in Tt.index:
    nneg, evmin = vq_stability(eng, K.loc[sn].values.astype(float), C.loc[sn].values.astype(float),
                               V.loc[sn].values.astype(complex), S.loc[sn].values.astype(int))
    rows.append(dict(snapshot=sn, n_bus_dVdQ_neg=nneg, eig_min_JR=evmin))
St = pd.DataFrame(rows).set_index("snapshot")
Tt = Tt.join(St); Tt["stable"] = (Tt.n_bus_dVdQ_neg == 0) & (Tt.eig_min_JR > 0)
Tt.to_csv(os.path.join(OUT, "v2_week_summary.csv"))
vm = V.map(abs); va = V.map(lambda z: np.degrees(np.angle(z)))
vm.to_csv(os.path.join(OUT, "v2_voltages_pu.csv")); va.to_csv(os.path.join(OUT, "v2_angles_deg.csv"))
Q.to_csv(os.path.join(OUT, "v2_regulators_Q_MVAr.csv")); S.to_csv(os.path.join(OUT, "v2_regulators_state.csv"))
K.to_csv(os.path.join(OUT, "v2_reactors_k.csv")); C.to_csv(os.path.join(OUT, "v2_capacitors_c.csv"))
ok = Tt.converged.astype(bool)
log(f"v2: Newton {int(Tt.newton_ok.sum())}/168, Q-limits {int(ok.sum())}, 0.90-1.10 {int((Tt.in_hard.astype(bool) & ok).sum())}, "
    f"0.95-1.05 {int((Tt.in_band.astype(bool) & ok).sum())}, stable {int(Tt.stable.sum())}, V {Tt.v_min.min():.3f}-{Tt.v_max.max():.3f}")

# v1 stability (stored solutions)
eng1n = R.load(os.path.join(HERE, "..", "reseau_315kv", "resultats", "reseau_315kv_final.nc"))
eng1 = R.ACEngine(eng1n, [g for g in eng1n.generators.index if g.startswith("REG ")])
rd = lambda f: pd.read_csv(os.path.join(HERE, "..", "reseau_315kv", "resultats", f), index_col=0, parse_dates=True)
V1, A1, S1 = rd("tensions_pu.csv"), rd("angles_deg.csv"), rd("regulateurs_etat.csv")
K1 = rd("inductances_k.csv").reindex(columns=eng1.buses).fillna(0); C1 = rd("condensateurs_c.csv").reindex(columns=eng1.buses).fillna(0)
rows = []
for sn in V1.index:
    Vc = V1.loc[sn, eng1.buses].values * np.exp(1j * np.radians(A1.loc[sn, eng1.buses].values))
    nneg, evmin = vq_stability(eng1, K1.loc[sn].values, C1.loc[sn].values, Vc, S1.loc[sn, eng1.regs].values.astype(int))
    rows.append(dict(snapshot=sn, n_bus_dVdQ_neg=nneg, eig_min_JR=evmin))
T1 = rd("semaine_resume.csv").join(pd.DataFrame(rows).set_index("snapshot"))
T1["stable"] = (T1.n_bus_dVdQ_neg == 0) & (T1.eig_min_JR > 0)
T1.to_csv(os.path.join(OUT, "v1_week_summary.csv"))
log(f"v1: stable {int(T1.stable.sum())}/168 ; unstable hours: {list(T1.index[~T1.stable].astype(str))}")

# ------------------------------------------------------------------ P-V curves at the peak (v1 and v2), shunts re-dispatched
def pv_curve(eng, label, lam_max=1.30, step=0.01, lam_min=0.90):
    P0, Q0 = eng.P.loc[PEAK].copy(), eng.Q.loc[PEAK].copy(); rows = []
    load = eng.n.loads_t.p_set.sum(axis=1)[PEAK]
    try:
        for lam in np.round(np.arange(lam_min, lam_max + 1e-9, step), 3):
            eng.P.loc[PEAK] = P0 * lam; eng.Q.loc[PEAK] = Q0 * lam
            r = R.solve_hour_reserve(eng, PEAK, max(0.0, min(R.K_LIGHT, R.K_PEAK - (lam - 1))), 1.0)
            if not r["ok"] or r["n_viol"] > 0:
                rows.append(dict(lam=lam, GW=load * lam / 1e3, solved=False)); break
            nneg, evmin = vq_stability(eng, r["k"], r["c"], r["V"], r["state"]); v = np.abs(r["V"])
            rows.append(dict(lam=lam, GW=load * lam / 1e3, solved=True, v_min=v.min(), bus_v_min=eng.buses[v.argmin()],
                             v_800=v[eng.ix["800"]] if "800" in eng.ix else np.nan, n_pq=int((r["state"] != 0).sum()),
                             stable=bool(nneg == 0 and evmin > 0), n_bus_dVdQ_neg=nneg, p_slack=r["p_slack"]))
            if lam > 1.0 and not rows[-1]["stable"]:
                break
    finally:
        eng.P.loc[PEAK], eng.Q.loc[PEAK] = P0, Q0
    D = pd.DataFrame(rows); D["case"] = label; return D

PV = pd.concat([pv_curve(eng1, "v1"), pv_curve(eng, "v2")])
PV.to_csv(os.path.join(OUT, "pv_curves_peak.csv"), index=False)
for c, d in PV.groupby("case"):
    st = d[d.solved & d.stable.fillna(False).astype(bool) & (d.lam >= 1.0)]
    log(f"{c}: last stable lambda at the peak = {st.lam.max() if len(st) else 'none'}")

# ------------------------------------------------------------------ validation of the v2 peak (PyPSA + pandapower) and MATPOWER export
W = dict(T=Tt, V=V, K=K, C=C, S=S)
D = R.validate_with_pypsa(n, eng, W, [PEAK, pd.Timestamp("2022-01-01 03:00")])
D.to_csv(os.path.join(OUT, "v2_validation_pypsa.csv")); log("PyPSA validation:\n" + D.round(8).to_string())
ppc = R.to_ppc(n, eng, W, PEAK)
try:
    pp = R.check_pandapower(ppc)
    vpp = pp.pop("v", None)
    pp["max_dV_vs_engine"] = float(np.abs(vpp - np.abs(V.loc[PEAK, eng.buses].values.astype(complex))).max()) if vpp is not None else np.nan
    pd.Series(pp).to_csv(os.path.join(OUT, "v2_validation_pandapower_peak.csv")); log("pandapower:", pp)
except Exception as e:
    log("pandapower failed:", repr(e))
R.write_matpower(ppc, os.path.join(OUT, "quebec_315kv_v2_peak.m"), name="quebec_315kv_v2_peak",
                 header="Quebec 315 kV+ network, 2022-01-03 17:00 (34.96 GW peak), v2 = run-of-river priority; stable V-Q solution. ECSE 458 capstone, Nov. 17 2026.")
with R.quiet():
    n.export_to_netcdf(os.path.join(OUT, "reseau_315kv_v2.nc"))

# ------------------------------------------------------------------ dispatch at the peak, v1 vs v2, by zone and carrier
def dispatch(nn):
    B = nn.buses; rows = []
    gp = nn.generators_t.p_set.loc[PEAK].reindex(nn.generators.index).fillna(0.0)
    for nm, r in nn.generators.iterrows():
        if r.carrier in ("load_shedding", "regulator") or "placeholder" in nm: continue
        rows.append(dict(name=nm, carrier=r.carrier, lat=B.y.get(r.bus), p_nom=r.p_nom, p=gp[nm]))
    sp = nn.storage_units_t.p_set.loc[PEAK].reindex(nn.storage_units.index).fillna(0.0)
    for nm, r in nn.storage_units.iterrows():
        rows.append(dict(name=nm, carrier="hydro reservoir", lat=B.y.get(r.bus), p_nom=r.p_nom, p=sp[nm]))
    d = pd.DataFrame(rows); d["zone"] = np.where(d.lat < 47.6, "south (< 47.6°N)", "north"); return d
d1 = dispatch(eng1n).assign(case="v1"); d2 = dispatch(n).assign(case="v2")
pd.concat([d1, d2]).to_csv(os.path.join(OUT, "dispatch_peak_v1_v2.csv"), index=False)
g = pd.concat([d1, d2]).groupby(["case", "zone", "carrier"]).agg(installed=("p_nom", "sum"), produced=("p", "sum")).round(0)
g.to_csv(os.path.join(OUT, "dispatch_peak_summary.csv")); log("\n" + g.to_string())

# ------------------------------------------------------------------ bus table for maps
B = n.buses.loc[eng.buses, ["x", "y", "v_nom"]].copy(); B.to_csv(os.path.join(OUT, "buses.csv"))
L = n.lines[["bus0", "bus1", "v_nom", "length", "num_parallel", "x", "b", "s_nom"]].copy(); L.to_csv(os.path.join(OUT, "lines.csv"))
log("done")
