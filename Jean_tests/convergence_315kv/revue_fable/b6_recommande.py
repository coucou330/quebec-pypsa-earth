"""B6 : la configuration recommandée = cas final + UNE correction : fil de l'eau prioritaire dans le LOPF
(coût marginal -0,2 $/MWh pour les centrales « ror »). On refait tout et on vérifie ce que la revue a appris
à vérifier : limites Q, bandes de tension, stabilité V-Q de chaque heure, marge stable à la pointe,
réactif fourni par des centrales à l'arrêt, puissance du bus slack."""
from b_outils import *
from b_variantes import q_engage, T0, _lopf0
import time
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _lopf0(n, solver=solver, log=log)
R.lopf = lopf_ror
with R.quiet(): n, info = R.build(T=T0)
R.lopf = _lopf0
eng = R.ACEngine(n, info["ac"]["regulators"])
t = time.time(); Tt, V, Q, K, C, S = R.solve_week_reserve(eng); print(f"semaine : {time.time() - t:.0f} s", flush=True)
ok = Tt.converged.astype(bool); load = n.loads_t.p_set.sum(axis=1)
def stable(Vc, st, k, c):
    Y = Ymat(eng, k, c)
    pv = sorted({eng.reg_bus[j] for j in range(len(eng.regs)) if st[j] == 0}); pvs = set(pv)
    pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]; pvpq = np.r_[pv, pq].astype(int)
    dth, dvm = eng._dSdV(Y, Vc)
    J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real], [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
    Ji = np.linalg.inv(J); m = len(pvpq); d = np.diag(Ji[m:, m:])
    ev = np.linalg.eigvals(np.linalg.inv(Ji[m:, m:])).real
    return int((d < 0).sum()), float(ev.min())
rows = []
for sn in Tt.index:
    nneg, evmin = stable(V.loc[sn].values.astype(complex), S.loc[sn].values.astype(int), K.loc[sn].values.astype(float), C.loc[sn].values.astype(float))
    rows.append(dict(heure=sn, bus_dVdQ_negatif=nneg, vp_min=evmin))
St = pd.DataFrame(rows).set_index("heure"); Tt["stable_VQ"] = (St.bus_dVdQ_negatif == 0) & (St.vp_min > 0)
Tt.to_csv("b6_semaine_resume.csv")
hard = Tt.in_hard.astype(bool) & ok; band = Tt.in_band.astype(bool) & ok
print(f"Newton {int(Tt.newton_ok.sum())}/168 | limites Q {int(ok.sum())}/168 | 0,90-1,10 : {int(hard.sum())} | 0,95-1,05 : {int(band.sum())} | stables V-Q : {int(Tt.stable_VQ.sum())}/168 | V {Tt.v_min.min():.3f}-{Tt.v_max.max():.3f}")
print(f"délestage LOPF {info['lopf']['shed_pct']:.4f} % | inductances enclenchées {Tt.q_react_MVAr.min() / 1e3:.1f} à {Tt.q_react_MVAr.max() / 1e3:.1f} GVAr | pertes (P slack) {Tt.p_slack.min():.0f} à {Tt.p_slack.max():.0f} MW")
cmp = sum(1 for sn in S.index for j, g in enumerate(eng.regs) if (S.at[sn, g] > 0 and abs(V.at[sn, eng.buses[eng.reg_bus[j]]]) > eng.vset[eng.reg_bus[j]] + 1e-3))
print("violations de complémentarité (> 0,001 pu) sur la semaine :", cmp)
# slack
su_b = n.storage_units[n.storage_units.bus == eng.slack_bus]; su_p = n.storage_units_t.p_set.reindex(columns=su_b.index).fillna(0).sum(axis=1)
tot = Tt.p_slack + su_p.reindex(Tt.index); print(f"bus slack : production max {tot.max():.0f} MW pour {su_b.p_nom.sum():.0f} MW installés ; heures au-dessus : {int((tot > su_b.p_nom.sum() + 1).sum())}/168")
# réactif de centrales à l'arrêt
viol = []
for sn in Tt.index:
    qmin, qmax, inst, on = q_engage(eng, sn); used = Q.loc[sn].values.astype(float)
    viol.append(dict(heure=sn, regs_au_dessus_du_Q_engage=int((used > qmax + 1).sum()), MVAr_au_dessus=float(np.clip(used - qmax, 0, None).sum())))
Vi = pd.DataFrame(viol).set_index("heure"); Vi.to_csv("b6_q_engage.csv")
print(f"réactif demandé à des machines à l'arrêt : {int((Vi.regs_au_dessus_du_Q_engage > 0).sum())}/168 heures ; à la pointe {Vi.MVAr_au_dessus[PEAK]:.0f} MVAr sur {Vi.regs_au_dessus_du_Q_engage[PEAK]} régulateurs ; max de la semaine {Vi.MVAr_au_dessus.max():.0f} MVAr")
# pointe avec limites Q « machines en marche »
qmin, qmax, inst, on = q_engage(eng, PEAK); q0 = (eng.qmin.copy(), eng.qmax.copy()); eng.qmin[:], eng.qmax[:] = np.minimum(qmin, 0), np.maximum(qmax, 0)
r = R.solve_hour_reserve(eng, PEAK, R.K_PEAK, 1.0); 
if r["ok"]:
    nneg, evmin = stable(r["V"], r["state"], r["k"], r["c"]); vm = np.abs(r["V"])
    print(f"pointe avec le réactif des seules machines en marche : Newton ok, régulateurs hors limites {r['n_viol']}, V {vm.min():.3f}-{vm.max():.3f}, stable V-Q : {nneg == 0 and evmin > 0}")
else: print("pointe avec le réactif des seules machines en marche : pas de solution")
# marge stable à la pointe, avec les deux jeux de limites Q
def marge(label):
    P0, Q0 = eng.P.loc[PEAK].copy(), eng.Q.loc[PEAK].copy(); rows = []
    for lam in np.round(np.arange(1.0, 1.301, 0.01), 3):
        eng.P.loc[PEAK] = P0 * lam; eng.Q.loc[PEAK] = Q0 * lam
        r = R.solve_hour_reserve(eng, PEAK, max(0.0, R.K_PEAK - (lam - 1)), 1.0)
        if not r["ok"] or r["n_viol"] > 0: break
        nneg, evmin = stable(r["V"], r["state"], r["k"], r["c"]); vm = np.abs(r["V"])
        rows.append(dict(lam=lam, GW=load[PEAK] * lam / 1e3, v_min=vm.min(), bus=eng.buses[vm.argmin()], n_pq=int((r["state"] != 0).sum()), stable=bool(nneg == 0 and evmin > 0)))
        if not rows[-1]["stable"]: break
    eng.P.loc[PEAK], eng.Q.loc[PEAK] = P0, Q0
    D = pd.DataFrame(rows); D.to_csv(f"b6_courbe_pv_{label}.csv", index=False)
    st = D[D.stable]; st90 = st[st.v_min >= 0.90]
    print(f"marge à la pointe ({label}) : dernier pas stable dans les limites Q : lambda {st.lam.max() if len(st) else float('nan')} ({st.GW.max() if len(st) else float('nan'):.1f} GW) ; avec V >= 0,90 : lambda {st90.lam.max() if len(st90) else float('nan')}")
    print(D.round(3).to_string(index=False))
marge("Q_machines_en_marche")
eng.qmin[:], eng.qmax[:] = q0
marge("Q_installe")
gs = n.generators[(n.generators.carrier == "ror")]; sud = gs[n.buses.y.reindex(gs.bus).values < 47.6]
print(f"fil de l'eau du sud à la pointe : {n.generators_t.p_set.loc[PEAK].reindex(sud.index).fillna(0).sum():.0f} MW produits pour {sud.p_nom.sum():.0f} MW installés")
