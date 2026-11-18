"""B3 : capacité réactive des centrales. Le modèle la calcule sur la puissance installée (p_nom), même pour
les machines à l'arrêt. Est-ce trop généreux ? On recalcule les limites Q de l'heure avec les seules machines
en marche (P > 1 MW dans le dispatch du LOPF), puis on refait l'heure."""
from b_outils import *
from b_variantes import q_engage
n, eng = final()
load = n.loads_t.p_set.sum(axis=1); lat = n.buses.y.reindex(eng.buses).values
sud = np.array([lat[b] < 47.6 for b in eng.reg_bus])
qmin0, qmax0 = eng.qmin.copy(), eng.qmax.copy()
Tt = rd("semaine_resume.csv"); Qr = rd("regulateurs_Q_MVAr.csv")
rows = []
for sn in [PEAK, pd.Timestamp("2022-01-03 16:00"), pd.Timestamp("2022-01-04 08:00"), pd.Timestamp("2022-01-02 16:00"), pd.Timestamp("2022-01-01 03:00")]:
    qmin, qmax, inst, on = q_engage(eng, sn)
    lamh = (load[sn] - load.min()) / (load.max() - load.min()); k0 = R.K_LIGHT - (R.K_LIGHT - R.K_PEAK) * lamh
    used = Qr.loc[sn, eng.regs].values
    row = dict(heure=sn, GW=round(load[sn] / 1e3, 1), machines_en_marche_pct=round(100 * on / inst), Qmax_installe=round(qmax0.sum()), Qmax_engage=round(qmax.sum()),
               Qmax_sud_installe=round(qmax0[sud].sum()), Qmax_sud_engage=round(qmax[sud].sum()), Q_fourni_sud=round(used[sud].sum()),
               regs_dont_Q_depasse_engage=int((used > qmax + 1).sum()))
    eng.qmin[:], eng.qmax[:] = np.minimum(qmin, 0), np.maximum(qmax, 0)
    r = R.solve_hour_reserve(eng, sn, k0, float(lamh))
    eng.qmin[:], eng.qmax[:] = qmin0, qmax0
    if r["ok"]:
        vm = np.abs(r["V"]); row.update(newton=True, n_viol=r["n_viol"], v_min=round(vm.min(), 3), v_max=round(vm.max(), 3), v_min_ref=round(Tt.v_min[sn], 3))
    else: row.update(newton=False)
    rows.append(row)
D = pd.DataFrame(rows); D.to_csv("b3_q_centrales.csv", index=False); print(D.T.to_string())
# qui fournit le Q du sud à la pointe ? centrales à l'arrêt ?
qmin, qmax, inst, on = q_engage(eng, PEAK); used = Qr.loc[PEAK, eng.regs].values
t = pd.DataFrame(dict(reg=eng.regs, lat=[round(lat[b], 1) for b in eng.reg_bus], Q=used.round(), Qmax_installe=qmax0.round(), Qmax_engage=qmax.round()))
t = t[(t.lat < 47.6)].sort_values("Qmax_installe", ascending=False); print("\nrégulateurs du sud à la pointe :"); print(t.head(20).to_string(index=False))
