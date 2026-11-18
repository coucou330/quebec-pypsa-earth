"""B4 (début) : que fait le LOPF des centrales du sud à la pointe ?"""
from b_outils import *
n, eng = final()
B = n.buses; g = n.generators; su = n.storage_units
gp = n.generators_t.p_set.loc[PEAK].reindex(g.index).fillna(0.0); sp = n.storage_units_t.p_set.loc[PEAK].reindex(su.index).fillna(0.0)
pmax = n.generators_t.p_max_pu.loc[PEAK].reindex(g.index).fillna(pd.Series(g.p_max_pu)) if len(n.generators_t.p_max_pu.columns) else g.p_max_pu
rows = []
for nm, r in g.iterrows():
    if r.carrier in ("load_shedding", "regulator") or "placeholder" in nm: continue
    rows.append(dict(nom=nm, bus=r.bus, carrier=r.carrier, lat=B.y.get(r.bus, np.nan), p_nom=r.p_nom, p=gp[nm], dispo=r.p_nom * float(pmax[nm]), cout=r.marginal_cost))
for nm, r in su.iterrows():
    rows.append(dict(nom=nm, bus=r.bus, carrier="hydro (réservoir)", lat=B.y.get(r.bus, np.nan), p_nom=r.p_nom, p=sp[nm], dispo=r.p_nom, cout=r.marginal_cost))
D = pd.DataFrame(rows); D["zone"] = np.where(D.lat < 47.6, "sud (< 47,6° N)", "nord"); D.to_csv("b4_dispatch_pointe.csv", index=False)
print(D.groupby(["zone", "carrier"]).agg(n=("p", "size"), installe=("p_nom", "sum"), disponible=("dispo", "sum"), produit=("p", "sum"), cout_min=("cout", "min"), cout_max=("cout", "max")).round(1).to_string())
ld = n.loads_t.p_set.loc[PEAK]; lz = ld.groupby(np.where(B.y.reindex(n.loads.bus.reindex(ld.index).values).values < 47.6, "sud", "nord")).sum()
print("\ncharge par zone (MW) :", lz.round().to_dict())
print("\ncentrales du sud de plus de 100 MW :"); print(D[(D.lat < 47.6) & (D.p_nom > 100)].sort_values("p_nom", ascending=False).round(1).to_string(index=False))
# sur la semaine : production des centrales du sud / disponible
sud_g = D[(D.lat < 47.6)].nom; gs = [x for x in sud_g if x in g.index]; ss = [x for x in sud_g if x in su.index]
ps = n.generators_t.p_set.reindex(columns=gs).fillna(0).sum(axis=1) + n.storage_units_t.p_set.reindex(columns=ss).fillna(0).sum(axis=1)
print(f"\nproduction du sud sur la semaine : min {ps.min():.0f}, moyenne {ps.mean():.0f}, max {ps.max():.0f} MW ; à la pointe {ps[PEAK]:.0f} MW pour {D[D.lat < 47.6].p_nom.sum():.0f} MW installés")
