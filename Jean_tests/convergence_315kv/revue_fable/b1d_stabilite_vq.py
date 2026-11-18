"""B1 (suite) : critère de stabilité V-Q sur les 168 solutions enregistrées dans ../reseau_315kv/resultats/.
Pour chaque bus PQ : dV_i/dQ_i (injecter du réactif au bus i doit faire MONTER sa tension). Un dV/dQ négatif
= la solution est sur la branche basse (instable) de la courbe : ce n'est pas un point de fonctionnement."""
from b_outils import *
n, eng = final()
V = rd("tensions_pu.csv"); A = rd("angles_deg.csv"); S = rd("regulateurs_etat.csv")
K = rd("inductances_k.csv").reindex(columns=eng.buses).fillna(0); C = rd("condensateurs_c.csv").reindex(columns=eng.buses).fillna(0)
load = n.loads_t.p_set.sum(axis=1); rows = []
for sn in V.index:
    Y = Ymat(eng, K.loc[sn].values, C.loc[sn].values)
    Vc = V.loc[sn, eng.buses].values * np.exp(1j * np.radians(A.loc[sn, eng.buses].values)); st = S.loc[sn, eng.regs].values.astype(int)
    pv = sorted({eng.reg_bus[j] for j in range(len(eng.regs)) if st[j] == 0}); pvs = set(pv)
    pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]; pvpq = np.r_[pv, pq].astype(int)
    dth, dvm = eng._dSdV(Y, Vc)
    J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real], [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
    Ji = np.linalg.inv(J); m = len(pvpq)
    d = np.diag(Ji[m:, m:]) * 100          # pu par 100 MVAr, pour chaque bus PQ
    neg = [(eng.buses[pq[i]], round(d[i], 4)) for i in np.where(d < 0)[0]]
    # valeurs propres de la matrice réduite J_R = dQ/dV (analyse modale V-Q de Gao-Morison-Kundur) : toutes > 0 si stable
    JR = np.linalg.inv(Ji[m:, m:]); ev = np.linalg.eigvals(JR).real
    rows.append(dict(heure=sn, GW=load[sn] / 1e3, v_min=V.loc[sn].min(), n_bus_dVdQ_negatif=len(neg), pire_dVdQ=d.max(), vp_min_JR=ev.min(), n_vp_negatives=int((ev < 0).sum()), bus_negatifs=str(neg[:6])))
D = pd.DataFrame(rows).set_index("heure"); D.to_csv("b1d_stabilite_vq.csv")
print("heures avec au moins un bus à dV/dQ négatif :", int((D.n_bus_dVdQ_negatif > 0).sum()), "/ 168 ; avec une valeur propre négative de J_R :", int((D.n_vp_negatives > 0).sum()))
print(D[(D.n_bus_dVdQ_negatif > 0) | (D.n_vp_negatives > 0)].round(4).to_string())
print("\nles 6 heures les plus chargées :"); print(D.sort_values("GW", ascending=False).head(6).drop(columns="bus_negatifs").round(4).to_string())
