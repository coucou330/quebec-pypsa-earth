"""C6 : résultat du réseau complet élagué (36 bus dont les charges restent au 315 kV, cf. c5) :
pointe et quelques heures avec limites Q et shunts manœuvrés, profil de tension par niveau, stabilité V-Q,
marge de charge. Comparaison avec le modèle 315 kV de même dispatch (fil de l'eau prioritaire)."""
from complet import *
import time, pickle
which = sys.argv[1] if len(sys.argv) > 1 else "ror"
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
excl = pickle.load(open("_c5_exclus_ror.pkl", "rb"))
n, info = build_complet(full, T, lopf=lopf_ror if which == "ror" else None, exclure=set(excl))
print(f"[{which}]", {k: v for k, v in info["complet"].items() if k != "bus_surcharges"}, flush=True)
eng = SparseEngine(n, info["ac"]["regulators"]); load = n.loads_t.p_set.sum(axis=1); nbrs = R.neighbours(eng)
print("bus :", eng.nb, "| condensateurs installés %.1f GVAr sur %d bus | inductances %.1f GVAr" % (eng.b_cap.sum() / 1e3, (eng.b_cap > 0).sum(), -eng.b_react.sum() / 1e3))
def stable(r):
    Y = sp.csr_matrix(eng.Y0 + np.diag(1j * (r["k"] * eng.b_react + r["c"] * eng.b_cap))); V = r["V"]; st = r["state"]
    pv = sorted({eng.reg_bus[j] for j in range(len(eng.regs)) if st[j] == 0}); pvs = set(pv)
    pq = np.array([i for i in range(eng.nb) if i != eng.sl and i not in pvs]); pvpq = np.r_[pv, pq].astype(int)
    I = Y @ V; dV = sp.diags(V); dI = sp.diags(I); dVn = sp.diags(V / np.abs(V))
    dth = (1j * dV @ (dI - Y @ dV).conj()).toarray(); dvm = (dV @ (Y @ dVn).conj() + dI.conj() @ dVn).toarray()
    J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real], [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
    Ji = np.linalg.inv(J); m = len(pvpq); d = np.diag(Ji[m:, m:])
    return int((d < 0).sum())
rows = []; HEURES = [PEAK, pd.Timestamp("2022-01-03 18:00"), pd.Timestamp("2022-01-04 10:00"), pd.Timestamp("2022-01-02 16:00"), pd.Timestamp("2022-01-05 04:00"), pd.Timestamp("2022-01-01 03:00")]
for sn in HEURES:
    lam = (load[sn] - load.min()) / (load.max() - load.min()); t = time.time()
    r = R.solve_hour_reserve(eng, sn, R.K_LIGHT - (R.K_LIGHT - R.K_PEAK) * lam, float(np.clip(lam, 0, 1)), nbrs=nbrs)
    if not r["ok"]: rows.append(dict(heure=sn, GW=load[sn] / 1e3, newton=False)); print(sn, "pas de solution", flush=True); continue
    tab, vm = bilan(n, eng, r); hv = vm[n.buses.v_nom.reindex(eng.buses).values >= 315]
    rows.append(dict(heure=sn, GW=round(load[sn] / 1e3, 1), newton=True, regs_hors_limites=r["n_viol"], regs_en_butee=int((r["state"] != 0).sum()), v_min=round(vm.min(), 3), v_max=round(vm.max(), 3),
                     bus_hors_090_110=int(((vm < 0.9) | (vm > 1.1)).sum()), bus_hors_095_105=int(((vm < 0.95) | (vm > 1.05)).sum()), v_min_315_735=round(hv.min(), 3), v_max_315_735=round(hv.max(), 3),
                     bus_dVdQ_negatif=stable(r), pertes_MW=round(r["p_slack"]), s=round(time.time() - t)))
    print(rows[-1], flush=True)
    if sn == PEAK:
        print(tab.to_string()); vm.to_csv(f"c6_tensions_pointe_{which}.csv")
        out = vm[(vm < 0.9) | (vm > 1.1)].sort_values(); print("bus hors 0,90-1,10 :", [(b, int(n.buses.v_nom[b]), round(v, 3)) for b, v in out.items()][:30])
pd.DataFrame(rows).to_csv(f"c6_heures_{which}.csv", index=False)
# marge à la pointe (Newton + limites Q + stabilité V-Q), pas de 1 %
P0, Q0 = eng.P.loc[PEAK].copy(), eng.Q.loc[PEAK].copy(); mr = []
for lam in np.round(np.arange(1.0, 1.31, 0.01), 3):
    eng.P.loc[PEAK] = P0 * lam; eng.Q.loc[PEAK] = Q0 * lam
    r = R.solve_hour_reserve(eng, PEAK, max(0.0, R.K_PEAK - (lam - 1)), 1.0, nbrs=nbrs)
    if not r["ok"] or r["n_viol"] > 0: break
    vm = np.abs(r["V"]); hv = vm[eng.vnom >= 315]; nneg = stable(r)
    mr.append(dict(lam=lam, GW=load[PEAK] * lam / 1e3, v_min=vm.min(), v_min_315_735=hv.min(), bus_sous_090=int((vm < 0.9).sum()), regs_en_butee=int((r["state"] != 0).sum()), stable=nneg == 0))
    print(mr[-1], flush=True)
    if nneg: break
eng.P.loc[PEAK], eng.Q.loc[PEAK] = P0, Q0
M = pd.DataFrame(mr); M.to_csv(f"c6_marge_{which}.csv", index=False)
if len(M): print(f"marge stable à la pointe : lambda {M.lam[M.stable].max()} ({M.GW[M.stable].max():.1f} GW)")
