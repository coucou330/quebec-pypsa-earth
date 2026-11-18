"""B1 : la marge nulle à la pointe est-elle un artefact de l'algorithme PV/PQ ?"""
from b_outils import *
import time
n, eng = final()
K = rd("inductances_k.csv").reindex(columns=eng.buses).fillna(0).loc[PEAK].values
C = rd("condensateurs_c.csv").reindex(columns=eng.buses).fillna(0).loc[PEAK].values
load = n.loads_t.p_set.sum(axis=1)[PEAK]
print(f"pointe : charge {load:.0f} MW ; shunts de la semaine : k moyen pondéré {(K * eng.b_react).sum() / eng.b_react.sum():.2f}, c {(C * eng.b_cap).sum() / eng.b_cap.sum():.2f}")
# 1. rappel : R.load_margin (départ à froid à chaque pas, critère de complémentarité 0,01)
for tol in (0.01, 9.0):
    d = R.load_margin(eng, PEAK, K, C, mode="fige", compl_tol=tol)
    print(f"R.load_margin figé, tolérance de complémentarité {tol} : lambda max {d.lam.max() if len(d) else float('nan')}")
# 2. continuation (départ à chaud, état conservé, un régulateur à la fois), shunts figés
t = time.time(); d = continuation(eng, PEAK, K, C); d.to_csv("b1_continuation_fige.csv", index=False)
print(f"continuation, shunts de la pointe figés : lambda max {d.lam.max():.4f} ({load * d.lam.max() / 1e3:.2f} GW), V min au dernier pas {d.v_min.iloc[-1]:.3f} (bus {d.bus.iloc[-1]}), "
      f"{d.n_pq.iloc[-1]} en butée, pire écart de complémentarité sur le parcours {d.compl.max():+.4f} pu ({time.time() - t:.0f} s)")
print(d.iloc[::max(1, len(d) // 8)].round(4).to_string(index=False))
# 3. continuation avec tous les condensateurs et moins d'inductances (ce qu'un exploitant ferait)
for kk in (0.30, 0.20, 0.10, 0.0):
    k2 = np.minimum(K, kk); c2 = np.ones(eng.nb)
    d2 = continuation(eng, PEAK, k2, c2)
    if len(d2): print(f"continuation, c = 1 partout, k plafonné à {kk} : lambda max {d2.lam.max():.4f} ; à lambda = 1 : V {d2.v_min.iloc[0]:.3f} min ; au dernier pas V min {d2.v_min.iloc[-1]:.3f} (bus {d2.bus.iloc[-1]}), compl {d2.compl.max():+.4f}")
    else: print(f"k plafonné à {kk} : pas de solution à lambda = 1")
