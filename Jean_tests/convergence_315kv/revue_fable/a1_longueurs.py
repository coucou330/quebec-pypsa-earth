"""A1 : la règle des longueurs est-elle bien appliquée dans le réseau final ?"""
from _commun import *
n, eng = final()
rep = pd.read_csv(os.path.join(RES, "longueurs_corrigees.csv"), index_col=0)
print(rep.statut.value_counts().to_string())
L = n.lines[n.lines.length > 0].copy()
B = n.buses
L["gc"] = hav(B.x[L.bus0].values, B.y[L.bus0].values, B.x[L.bus1].values, B.y[L.bus1].values)
L["ratio"] = L.length / L.gc
L["statut"] = rep.statut.reindex(L.index)
L["methode"] = rep.methode.reindex(L.index) if "methode" in rep else ""
print("\nlongueur / vol d'oiseau, par statut :")
print(L.groupby("statut", dropna=False).ratio.describe()[["count", "min", "50%", "max"]].round(3).to_string())
bad = L[(L.ratio < 0.999) | (L.ratio > 1.7)]
print(f"\nlignes < vol d'oiseau ou > 1,7 : {len(bad)} sur {len(L)}")
cols = ["bus0", "bus1", "v_nom", "length", "gc", "ratio", "num_parallel", "statut"]
print(bad.sort_values("ratio")[cols].round(2).to_string())
# poids de ces lignes
for st, g in L.groupby("statut", dropna=False):
    print(f"{st}: {len(g)} lignes, {(g.length*g.num_parallel).sum():.0f} circuit-km, dont 735 kV {(g.length*g.num_parallel)[g.v_nom>=700].sum():.0f}")
# lignes non corrigées : facteur Mercator attendu à leur latitude
nc = L[L.statut.isin(["modifiée depuis elec_full", "sans géométrie"]) | L.statut.isna()].copy()
nc["sec_lat"] = 1 / np.cos(np.radians((B.y[nc.bus0].values + B.y[nc.bus1].values) / 2))
nc["ratio_sur_sec"] = nc.ratio / nc.sec_lat
print("\nnon corrigées : ratio longueur/vol d'oiseau comparé à 1/cos(lat) (≈1,0-1,1 => encore en Mercator)")
print(nc.sort_values("length", ascending=False)[["bus0","bus1","v_nom","length","gc","ratio","sec_lat","ratio_sur_sec","statut"]].round(2).head(40).to_string())
L.to_csv("a1_longueurs.csv")
