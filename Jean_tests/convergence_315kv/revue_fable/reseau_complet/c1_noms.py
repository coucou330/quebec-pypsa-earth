import sys, os; sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _commun import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))
print("charges solved (5):", list(src.loads.index[:5]), "| full (5):", list(full.loads.index[:5]))
common = src.loads.index.intersection(full.loads.index); print("noms communs :", len(common), "/", len(src.loads))
print(src.loads.head(3)[["bus", "p_set"]].to_string()); print(full.loads.head(3)[["bus", "p_set"]].to_string())
if len(common):
    ob = full.loads.bus[common]; print("tension du bus d'origine des charges communes :", full.buses.v_nom[ob.values].value_counts().sort_index().to_dict())
    moved = (ob.values != src.loads.bus[common].values).sum(); print("charges déplacées par la réduction :", moved)
print("gén. solved (8):", list(src.generators.index[:8])); print("gén. full (5):", list(full.generators.index[:5]), "| communs", len(src.generators.index.intersection(full.generators.index)), "| stockages communs", len(src.storage_units.index.intersection(full.storage_units.index)), "/", len(src.storage_units))
red = R.load(os.path.join(R.CURRENT, "elec_reduced.nc")); print("elec_reduced :", len(red.buses), "bus,", len(red.loads), "charges, noms communs avec full :", len(red.loads.index.intersection(full.loads.index)), "; avec solved :", len(red.loads.index.intersection(src.loads.index)))
print("gén. reduced communs avec full :", len(red.generators.index.intersection(full.generators.index)), "/", len(red.generators), "; stockages", len(red.storage_units.index.intersection(full.storage_units.index)), "/", len(red.storage_units))
