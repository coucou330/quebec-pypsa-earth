"""B1 (suite) : marge à la pointe avec des régulateurs à statisme (pas de bascule PV/PQ).
On part de 85 % de la charge de pointe et on monte par pas de 0,25 %, shunts de la pointe figés."""
from b_outils import *
n, eng = final()
K = rd("inductances_k.csv").reindex(columns=eng.buses).fillna(0).loc[PEAK].values
C = rd("condensateurs_c.csv").reindex(columns=eng.buses).fillna(0).loc[PEAK].values
load = n.loads_t.p_set.sum(axis=1)[PEAK]
for slope in (0.005, 0.01, 0.02, 0.04):
    d = marge_statisme(eng, PEAK, K, C, slope=slope, lam_min=0.85); d.to_csv(f"b1_pv_statisme_{slope}.csv", index=False)
    if not len(d): print(f"statisme {slope:.1%} : pas de solution à lambda = 0,85"); continue
    s = resume_marge(d, load)
    print(f"statisme {slope:.1%} : nez de la courbe à lambda {s['lam_nez']:.4f} ({s['GW_nez']} GW) ; lambda max avec V >= 0,90 : {s['lam_090']:.4f} ; "
          f"au nez : V min {d.v_min.iloc[-1]:.3f} (bus {d.bus.iloc[-1]}), {d.sat_hi.iloc[-1]} régulateurs à Qmax, V max {d.v_max.max():.3f} ; à lambda = 1 : V min {s['vmin_1']:.3f}")
d = pd.read_csv("b1_pv_statisme_0.02.csv"); print(d.iloc[::max(1, len(d) // 10)].round(4).to_string(index=False))
