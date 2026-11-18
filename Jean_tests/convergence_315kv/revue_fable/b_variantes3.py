"""B2 (suite) : Kron étendu au sud, à dispatch corrigé (fil de l'eau prioritaire, B4), pour séparer l'effet
de Kron de celui du dispatch. La boîte « Sud » (lat 44,9-47,6 ; lon -77,5 à -69,0) contient aussi des bus de
l'Ontario, de l'État de New York et du Maine : le réseau d'elec_full y est relié en AC, et Kron les
transforme en chemins parallèles qui n'existent pas en exploitation. Variante « Québec seulement » : ces bus
sont sortis de la boîte (frontière approchée par une ligne brisée : rivière des Outaouais, 45e parallèle, Maine).
Usage : python b_variantes3.py i   (0 : boîte complète ; 1 : Québec seulement ; 2 : Québec seulement, tout le Québec)"""
import sys, copy
import b_variantes as BV
from b_variantes import *

def hors_quebec(B):
    lon, lat = B.x.values, B.y.values
    west = np.interp(lon, [-80.0, -79.5, -78.7, -77.5, -77.1, -76.3, -75.7, -74.4], [47.0, 46.6, 46.3, 46.05, 45.83, 45.55, 45.42, 45.57])
    east = np.interp(lon, [-74.4, -71.5, -70.6, -70.0, -69.3, -69.0, -67.8], [45.0, 45.0, 45.3, 46.0, 47.3, 47.45, 47.1])
    border = np.where(lon < -74.4, west, east)
    return pd.Series((lat < border) & (lon < -67.8), index=B.index)

full_q = R.load(os.path.join(R.CURRENT, "elec_full.nc"))
hq = hors_quebec(full_q.buses); print("bus d'elec_full hors Québec (approx.) sous la frontière sud :", int(hq.sum()), flush=True)
full_q.buses.loc[hq, "y"] = -99.0          # hors de toute boîte de Kron
_kron0 = R.kron_equivalent
def kron_quebec(full, model, T, box=R.SAGUENAY_BOX, **kw):
    return _kron0(full_q, model, T, box=box, **kw)
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return BV._lopf0(n, solver=solver, log=log)
SUD = {"Saguenay": R.SAGUENAY_BOX, "Sud": dict(lat=(44.9, 47.6), lon=(-77.5, -69.0))}
QC = {"Québec": dict(lat=(44.0, 63.0), lon=(-80.0, -56.0))}
VAR3 = [
    ("B4 + B2 : fil de l'eau prioritaire, Kron Saguenay + sud (boîte complète, Ontario/É.-U. inclus)", dict(lopf=lopf_ror, regions=SUD), False),
    ("B4 + B2 : fil de l'eau prioritaire, Kron Saguenay + sud, Québec seulement", dict(lopf=lopf_ror, regions=SUD), True),
    ("B4 + B2 : fil de l'eau prioritaire, Kron tout le Québec, Québec seulement", dict(lopf=lopf_ror, regions=QC), True),
]
if __name__ == "__main__":
    i = int(sys.argv[1]); label, opt, qc = VAR3[i]
    src = R.load(os.path.join(R.CURRENT, "elec_solved.nc")); BV.eng_ref_load = src.loads_t.p_set.sum(axis=1).reindex(BV.REF.index)
    R.kron_equivalent = kron_quebec if qc else _kron0
    try: row = BV.run(label, opt, None)
    except Exception as e:
        import traceback; traceback.print_exc(); row = dict(variante=label, erreur=repr(e)[:300])
    print(row, flush=True)
    out = f"b_variantes_part3{i}.csv"; pd.DataFrame([row]).to_csv(out, index=False)
