"""Import de reseau315 (dossier voisin, non modifié) + chargement du réseau final."""
import os, sys, warnings
warnings.filterwarnings("ignore")
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "reseau_315kv"))
import numpy as np, pandas as pd
import reseau315 as R
RES = os.path.join(HERE, "..", "reseau_315kv", "resultats")
PEAK = pd.Timestamp("2022-01-03 17:00")
def final():
    n = R.load(os.path.join(RES, "reseau_315kv_final.nc"))
    regs = [g for g in n.generators.index if g.startswith("REG ")]
    return n, R.ACEngine(n, regs)
def rd(f):
    return pd.read_csv(os.path.join(RES, f), index_col=0, parse_dates=True)
def hav(lon1, lat1, lon2, lat2):
    lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
    a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
    return 2 * 6371.0 * np.arcsin(np.sqrt(a))
