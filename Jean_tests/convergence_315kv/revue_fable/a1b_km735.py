"""A1b : circuit-km à 735 kV dans le réseau final, comparés aux 12 581 km d'HQ (rapport annuel 2023)."""
from _commun import *
n, eng = final()
L = n.lines[(n.lines.v_nom >= 700) & (n.lines.length > 0)]
print("735 kV : lignes", len(L), "| somme longueurs", round(L.length.sum()), "km | circuit-km (× num_parallel)", round((L.length * L.num_parallel).sum()))
print(L.num_parallel.value_counts().to_string())
lab = L[L.bus0.isin(["469"]) | L.bus1.isin(["469"])]
print("dont Churchill Falls-Montagnais (Labrador/Québec) :", round((lab.length * lab.num_parallel).sum()), "km")
src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))
S = src.lines.reindex(L.index)
d = pd.DataFrame(dict(bus0=L.bus0, bus1=L.bus1, km=L.length.round(), np_final=L.num_parallel, np_source=S.num_parallel))
ch = d[d.np_final != d.np_source]
print(f"\nnum_parallel changé par E3 (script de Nathan) : {len(ch)} lignes 735 kV, +{((ch.np_final-ch.np_source)*ch.km).sum():.0f} circuit-km")
print(ch.sort_values("km", ascending=False).to_string())
# doublons : plusieurs lignes entre la même paire de bus
pair = L.apply(lambda r: tuple(sorted([r.bus0, r.bus1])), axis=1)
g = pd.DataFrame(dict(pair=pair, km=L.length, np=L.num_parallel)).groupby("pair").agg(n=("km", "size"), km=("km", "mean"), circuits=("np", "sum"))
print("\ncorridors (paires de bus) avec le plus de circuits :"); print(g.sort_values("circuits", ascending=False).head(15).round(1).to_string())
print("circuits par corridor :", g.circuits.value_counts().sort_index().to_dict())
