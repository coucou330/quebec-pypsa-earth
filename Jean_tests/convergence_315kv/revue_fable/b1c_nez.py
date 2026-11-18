"""B1 (suite) : la pointe est-elle au nez de la courbe P-V ? Deux indicateurs indépendants de la logique PV/PQ :
(a) courbe P-V de 85 % à 102 % de la charge de pointe (shunts re-manœuvrés à chaque pas) ;
(b) plus petite valeur singulière du jacobien de Newton à la solution (0 = nez), comparée à d'autres heures ;
(c) sensibilité dV/dQ au bus le plus bas (bus 800) : explose près du nez."""
from b_outils import *
n, eng = final()
load = n.loads_t.p_set.sum(axis=1)
def jac(eng, Y, V, state):
    pv = sorted({eng.reg_bus[j] for j in range(len(eng.regs)) if state[j] == 0}); pvs = set(pv)
    pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]; pvpq = np.r_[pv, pq].astype(int)
    dth, dvm = eng._dSdV(Y, V)
    J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real], [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
    return J, pvpq, pq
def indicateurs(r):
    Y = Ymat(eng, r["k"], r["c"]); J, pvpq, pq = jac(eng, Y, r["V"], r["state"])
    sv = np.linalg.svd(J, compute_uv=False)
    i800 = eng.ix["800"]; out = dict(sv_min=sv[-1], sv_rel=sv[-1] / sv[0])
    if i800 in pq:
        e = np.zeros(J.shape[0]); e[len(pvpq) + pq.index(i800)] = 1.0     # +1 MVAr injecté au bus 800
        dx = np.linalg.solve(J, e); out["dV800_dQ800_pu_par_100MVAr"] = 100 * dx[len(pvpq) + pq.index(i800)]
    return out
P0, Q0 = eng.P.loc[PEAK].copy(), eng.Q.loc[PEAK].copy(); rows = []
for lam in [0.85, 0.90, 0.94, 0.96, 0.98, 0.99, 0.995, 1.0, 1.005, 1.01, 1.015, 1.02]:
    eng.P.loc[PEAK] = P0 * lam; eng.Q.loc[PEAK] = Q0 * lam
    r = R.solve_hour_reserve(eng, PEAK, max(0.0, min(0.85, R.K_PEAK - (lam - 1) * 3)), 1.0)
    if not r["ok"]: rows.append(dict(lam=lam, ok=False)); continue
    vm = np.abs(r["V"]); cv = R.complementarity_violations(eng, r, 0.001)
    rows.append(dict(lam=lam, GW=load[PEAK] * lam / 1e3, ok=True, n_viol=r["n_viol"], compl=max([c[1] for c in cv] + [0]), v_min=vm.min(), bus=eng.buses[vm.argmin()], n_pq=int((r["state"] != 0).sum()), **indicateurs(r)))
eng.P.loc[PEAK], eng.Q.loc[PEAK] = P0, Q0
D = pd.DataFrame(rows); D.to_csv("b1c_courbe_pv_pointe.csv", index=False); print(D.round(5).to_string(index=False))
print("\nmêmes indicateurs à d'autres heures de la semaine (charge telle quelle) :")
rows = []
for sn in [pd.Timestamp("2022-01-01 03:00"), pd.Timestamp("2022-01-02 16:00"), pd.Timestamp("2022-01-03 09:00"), pd.Timestamp("2022-01-03 16:00"), PEAK, pd.Timestamp("2022-01-04 10:00"), pd.Timestamp("2022-01-06 18:00")]:
    lamh = (load[sn] - load.min()) / (load.max() - load.min())
    r = R.solve_hour_reserve(eng, sn, R.K_LIGHT - (R.K_LIGHT - R.K_PEAK) * lamh, float(lamh))
    vm = np.abs(r["V"]); rows.append(dict(heure=sn, GW=load[sn] / 1e3, v_min=vm.min(), n_pq=int((r["state"] != 0).sum()), **indicateurs(r)))
print(pd.DataFrame(rows).round(5).to_string(index=False))
