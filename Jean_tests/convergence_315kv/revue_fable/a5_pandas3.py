"""A5 : le module tourne-t-il avec pandas 3.x ? (PYTHONPATH=/tmp/pd3 : pandas 3 installé à part)"""
from _commun import *
import pandas, pypsa, numpy, time
print("pandas", pandas.__version__, "| numpy", numpy.__version__, "| pypsa", pypsa.__version__, flush=True)
t = time.time()
n, info = R.build()
print("build ok", round(time.time() - t), "s ; délestage LOPF %.4f %%" % info["lopf"]["shed_pct"], flush=True)
eng = R.ACEngine(n, info["ac"]["regulators"])
sns = list(n.snapshots[[0, 3, 40, 65, 66, 100, 167]])
Tt, V, Q, K, C, S = R.solve_week_reserve(eng, sns)
print(Tt[["converged", "regulateurs_hors_limites", "v_min", "v_max", "n_pq", "q_react_MVAr"]].round(4).to_string())
ref = rd("semaine_resume.csv").loc[sns]
print("écart avec resultats/ (pandas 2.3.3) : v_min %.1e, v_max %.1e, q_react %.1e MVAr" % ((Tt.v_min - ref.v_min).abs().max(), (Tt.v_max - ref.v_max).abs().max(), (Tt.q_react_MVAr - ref.q_react_MVAr).abs().max()))
W = dict(T=Tt, V=V, K=K, C=C, S=S)
D = R.validate_with_pypsa(n, eng, W, sns[:2]); print(D.round(8).to_string())
ppc = R.to_ppc(n, eng, W, sns[3]); print("to_ppc ok", ppc["bus"].shape, ppc["gen"].shape)
