"""C0 : inventaire d'elec_full.nc (réseau complet, toutes tensions) et lien avec elec_solved.nc (315 kV)."""
import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _commun import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))
print("snapshots full:", len(full.snapshots), full.snapshots[0], "| solved:", len(src.snapshots))
for nm, n in (("elec_full", full), ("elec_solved", src)):
    print(f"\n== {nm} : {len(n.buses)} bus, {len(n.lines)} lignes, {len(n.transformers)} transfos, {len(n.links)} liens, {len(n.generators)} gén., {len(n.storage_units)} stockages, {len(n.loads)} charges, {len(n.shunt_impedances)} shunts")
    print("bus par tension :", n.buses.v_nom.value_counts().sort_index().to_dict())
    print("lignes par tension :", n.lines.v_nom.value_counts().sort_index().to_dict() if "v_nom" in n.lines else "")
    print("carriers bus :", n.buses.carrier.value_counts().to_dict())
    print("charges : p_set horaire", n.loads_t.p_set.shape, "pointe %.0f MW" % (n.loads_t.p_set.sum(axis=1).max() if n.loads_t.p_set.size else n.loads.p_set.sum()))
    print("gén. carriers :", n.generators.groupby("carrier").p_nom.sum().round().to_dict())
    print("generators_t.p :", n.generators_t.p.shape, "p_set", n.generators_t.p_set.shape, "p_max_pu", n.generators_t.p_max_pu.shape)
full.determine_network_topology()
sn = full.buses.sub_network.value_counts(); print("\nîlots (sous-réseaux) d'elec_full :", len(sn), "; tailles :", sn.head(8).to_dict())
main = sn.index[0]; bm = full.buses.index[full.buses.sub_network == main]
print("îlot principal : bus par tension", full.buses.v_nom[bm].value_counts().sort_index().to_dict())
ld = full.loads; print("charges dans l'îlot principal : %d / %d" % (ld.bus.isin(bm).sum(), len(ld)))
print("bus d'elec_solved absents d'elec_full :", [b for b in src.buses.index if b not in full.buses.index][:10])
print("bus d'elec_solved hors îlot principal d'elec_full :", int((~src.buses.index[src.buses.index.isin(full.buses.index)].isin(bm)).sum()))
T = full.transformers; print("\ntransfos : s_nom", T.s_nom.describe()[["min", "50%", "max"]].round().to_dict(), "| x", T.x.describe()[["min", "50%", "max"]].round(3).to_dict(), "| r", T.r.describe()[["min", "max"]].to_dict(), "| types", T.type.value_counts().head(4).to_dict())
L = full.lines; print("lignes : types", L.type.value_counts().head(8).to_dict()); print("x nul :", int((L.x == 0).sum()), "| longueur nulle :", int((L.length == 0).sum()), "| num_parallel 0 :", int((L.num_parallel == 0).sum()))
print("liens :", full.links.carrier.value_counts().to_dict(), full.links[["bus0", "bus1", "p_nom"]].head(8).to_string())
# où sont les charges et la production d'elec_solved par rapport à elec_full ?
print("\ncharges full par tension de bus :", full.loads.assign(v=full.buses.v_nom[full.loads.bus].values).groupby("v").size().to_dict())
print("gén. full par tension de bus (MW) :", full.generators.assign(v=full.buses.v_nom[full.generators.bus].values).groupby("v").p_nom.sum().round().to_dict())
print("stockages full :", full.storage_units.assign(v=full.buses.v_nom[full.storage_units.bus].values).groupby("v").p_nom.sum().round().to_dict())
print("charges solved par tension :", src.loads.assign(v=src.buses.v_nom[src.loads.bus].values).groupby("v").size().to_dict())
