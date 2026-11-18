# -*- coding: utf-8 -*-
"""
reseau315.py — un réseau 315 kV du Québec qui converge en AC power flow, avec un réactif réaliste.

Capstone ECSE 458 (Jean Flichy et équipe). Construit par Claude (conversation claude.ai) les 16-17 nov.
2026 à partir du modèle PyPSA-Earth de Nathan Guan. Explications, physique, hypothèses et résultats :
notebook `Reseau_315kV_converge.ipynb` (même dossier). À relire et refaire étape par étape.

Point de départ : `network/networks_current/elec_solved.nc` (315 kV de Nathan, LOPF déjà configuré).
Rien n'est écrit dans `network/networks_current/`.

Étapes (une fonction chacune ; `build()` les enchaîne, `solve_week_reserve()` résout la semaine) :
  E1   fix_lengths           longueurs gonflées par la projection Web Mercator (x, b, r trop grands de ~50 %)
  E2   fix_churchill_falls   poste de Churchill Falls coupé en deux bus -> 6 circuits au lieu de 3
  E2b  kron_equivalent       réseau < 315 kV du Saguenay supprimé par la réduction -> équivalent de Kron
  E3   correct_circuits      nombre de circuits des corridors de départ (script de Nathan, inchangé)
  E4   series_compensation   40 % sur les corridors où HQ a vraiment de la compensation série
  E5   st_clair_ratings,     3 × St Clair avec les vraies longueurs + plafond de stabilité (V²/x) sin 40°,
       lopf                  puis LOPF relancé
  E6   prepare_ac            charges cos phi 0,95 ; un régulateur de tension par bus (centrales + SVC/CS
                             réels d'HQ, limites Q) ; inductances de ligne commutables ; consigne 1,02 pu
  E6c  capacitor_banks_from_load   bancs de condensateurs commutables aux postes de charge
  E7   ACEngine + solve_week_reserve   Newton-Raphson heure par heure, bascule PV -> PQ aux limites Q,
                             commutation des shunts pour garder de la réserve dans les régulateurs (logique
                             d'exploitation d'HQ), puis vérification par le solveur de PyPSA

Compatible PyPSA 1.3.0 (venv `pypsa-env`) ; testé sous Linux avec PyPSA 1.3.0, HiGHS 1.16, pandas 2.3/3.0.
"""
import contextlib
import io
import logging
import os
import sys
import time
import warnings

import numpy as np
import pandas as pd
import pypsa

warnings.filterwarnings("ignore")
logging.getLogger("pypsa").setLevel(logging.ERROR)
logging.getLogger("linopy").setLevel(logging.ERROR)


# ============================================================================ chemins
def find_repo(start=None):
    """Remonte depuis ce fichier jusqu'au dépôt (dossier qui contient network/run_pf.py)."""
    d = os.path.abspath(start or os.path.dirname(os.path.abspath(__file__)))
    for _ in range(8):
        if os.path.isfile(os.path.join(d, "network", "run_pf.py")):
            return d
        d = os.path.dirname(d)
    raise FileNotFoundError("Dépôt quebec-pypsa-earth introuvable (network/run_pf.py).")

REPO = os.environ.get("QC_REPO") or find_repo()
NETWORK_DIR = os.path.join(REPO, "network")
CURRENT = os.path.join(NETWORK_DIR, "networks_current")
if NETWORK_DIR not in sys.path:
    sys.path.insert(0, NETWORK_DIR)

import run_pf  # noqa: E402  (préparation AC de Nathan)
import correct_sending_end_circuits as _csc  # noqa: E402
from st_clair import load_line_params, st_clair_s_nom_mva  # noqa: E402
from apply_hq_line_characteristics import VOLTAGE_TIER_MAP, load_hq_params  # noqa: E402


@contextlib.contextmanager
def quiet():
    with contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(io.StringIO()):
        yield

def load(path):
    with quiet():
        return pypsa.Network(path)



# ============================================================================ E1 longueurs (Web Mercator)
def mercator_table(full):
    """Pour chaque ligne qui a une géométrie OSM : longueur géodésique (vraie) et longueur
    de la même géométrie en EPSG:3857 (ce que PyPSA-Earth a stocké)."""
    from pyproj import Transformer
    from shapely import wkt
    tr = Transformer.from_crs("EPSG:4326", "EPSG:3857", always_xy=True)

    def hav(lon1, lat1, lon2, lat2):
        lon1, lat1, lon2, lat2 = map(np.radians, [lon1, lat1, lon2, lat2])
        a = np.sin((lat2 - lat1) / 2) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2) ** 2
        return 2 * 6371.0 * np.arcsin(np.sqrt(a))

    rows = {}
    L = full.lines
    has = L.geometry.astype(str).str.startswith(("LINESTRING", "MULTILINESTRING"))
    for name in L.index[has]:
        g = wkt.loads(str(L.at[name, "geometry"]))
        parts = g.geoms if g.geom_type == "MultiLineString" else [g]
        t = m = 0.0
        lat_w = 0.0
        for p in parts:
            c = np.asarray(p.coords)
            seg = hav(c[:-1, 0], c[:-1, 1], c[1:, 0], c[1:, 1])
            t += seg.sum()
            lat_w += (seg * (c[:-1, 1] + c[1:, 1]) / 2).sum()
            x, y = tr.transform(c[:, 0], c[:, 1])
            m += np.hypot(np.diff(x), np.diff(y)).sum() / 1e3
        b0, b1 = L.at[name, "bus0"], L.at[name, "bus1"]
        x0, y0, x1, y1 = full.buses.x[b0], full.buses.y[b0], full.buses.x[b1], full.buses.y[b1]
        gc = hav(x0, y0, x1, y1)
        c0, c1 = parts[0].coords[0], parts[-1].coords[-1]
        mis = min(hav(c0[0], c0[1], x0, y0) + hav(c1[0], c1[1], x1, y1), hav(c0[0], c0[1], x1, y1) + hav(c1[0], c1[1], x0, y0))
        rows[name] = dict(v_nom=L.at[name, "v_nom"], length_stored=L.at[name, "length"], true_km=t,
                          merc_km=m, gc_km=gc, end_mismatch_km=mis, lat=lat_w / t if t else np.nan)
    T = pd.DataFrame.from_dict(rows, orient="index")
    T["factor"] = T.true_km / T.merc_km
    T["stored_over_true"] = T.length_stored / T.true_km
    T["stored_over_merc"] = T.length_stored / T.merc_km
    T["sec_lat"] = 1 / np.cos(np.radians(T.lat))
    T["sinuosity"] = T.true_km / T.gc_km
    T["reliable"] = (T.end_mismatch_km <= 5.0) & T.sinuosity.between(1.0, 1.6)
    return T

def fix_lengths(n, T, log=None):
    """Longueur vraie de chaque ligne (puis r, x, b à l'échelle, car tous proportionnels à la longueur) :
    - tracé OSM fiable (extrémités à <= 5 km des bus, sinuosité tracé / vol d'oiseau entre 1,0 et 1,6) :
      longueur géodésique du tracé ;
    - sinon : vol d'oiseau × sinuosité médiane des tracés fiables (le tracé ou la longueur stockée ne
      correspondent pas à la ligne : segment droit mesuré en Mercator, tracé qui fait un détour...).
    On saute les lignes sans géométrie, celles modifiées depuis elec_full (réduction, ajouts de Nathan) et
    celles déjà à moins de 2 % de la bonne longueur."""
    sin_med = float(T.loc[T.reliable, "sinuosity"].median())
    rep = []
    for name in n.lines.index:
        if name not in T.index:
            rep.append(dict(line=name, statut="sans géométrie")); continue
        L0 = n.lines.at[name, "length"]
        if abs(L0 / T.at[name, "length_stored"] - 1) >= 0.01:
            rep.append(dict(line=name, statut="modifiée depuis elec_full")); continue
        if T.at[name, "reliable"]:
            L1, how = T.at[name, "true_km"], "tracé géodésique"
        else:
            L1, how = sin_med * T.at[name, "gc_km"], "vol d'oiseau × sinuosité"
        if abs(L1 / L0 - 1) < 0.02:
            rep.append(dict(line=name, statut="déjà à la bonne longueur", methode=how)); continue
        f = L1 / L0
        for c in ["length", "r", "x", "b"]:
            n.lines.at[name, c] *= f
        rep.append(dict(line=name, statut="corrigée", methode=how, v_nom=n.lines.at[name, "v_nom"], km_avant=L0,
                        km_apres=L1, facteur=f, mercator=T.at[name, "factor"]))
    n.calculate_dependent_values()
    rep = pd.DataFrame(rep).set_index("line")
    rep.attrs["sinuosite_mediane"] = sin_med
    if log is not None:
        ok = rep[rep.statut == "corrigée"]
        log.append(("E1 longueurs", f"{len(ok)} lignes corrigées ({(ok.methode == 'tracé géodésique').sum()} par leur tracé, "
                    f"{(ok.methode != 'tracé géodésique').sum()} par vol d'oiseau × {sin_med:.2f}), facteur médian {ok.facteur.median():.3f}"))
    return rep

# ============================================================================ E2 Churchill Falls
CF_PLANT_BUS = "465-735kv"   # bus ajouté par add_churchill_falls_tie.py (centrale + 3 lignes L7051-L7053)

CF_OSM_BUS = "469"           # vrai poste 735 kV OSM de Churchill Falls (3 lignes 1568-1570 vers Montagnais)

CF_DUPLICATE_LINES = ["465-735kv-306 L7051", "465-735kv-306 L7052", "465-735kv-306 L7053"]

def fix_churchill_falls(n, log=None):
    """Les 3 lignes OSM 1568-1570 (469-306) et les 3 lignes ajoutées L7051-L7053 (465-735kv-306)
    sont les mêmes 3 circuits réels Churchill Falls-Montagnais. On garde les lignes OSM, on rattache
    la centrale (et la charge locale) au poste OSM 469 et on supprime le doublon."""
    if CF_PLANT_BUS not in n.buses.index:
        return 0
    for c in ["storage_units", "generators", "loads"]:
        df = getattr(n, c)
        idx = df.index[df.bus == CF_PLANT_BUS]
        df.loc[idx, "bus"] = CF_OSM_BUS
    dup = [l for l in CF_DUPLICATE_LINES if l in n.lines.index]
    with quiet():
        n.remove("Line", dup)
        n.remove("Bus", CF_PLANT_BUS)
    if log is not None:
        log.append(("E2 Churchill Falls", f"{len(dup)} lignes en double supprimées, centrale rattachée au bus {CF_OSM_BUS}"))
    return len(dup)



# ============================================================================ E2b équivalent de Kron régional (Saguenay)
SAGUENAY_BOX = dict(lat=(47.6, 49.0), lon=(-73.0, -70.3))
MONTREAL_BOX = dict(lat=(45.25, 45.85), lon=(-74.2, -73.25))
KRON_REGIONS = {"Saguenay": SAGUENAY_BOX, "Montréal": MONTREAL_BOX}

BETA_50HZ = 1.2   # types de lignes PyPSA-Earth < 315 kV : x et b calculés à 50 Hz (sous-estimés de 17 %, cf. cours)

def kron_equivalent(full, model, T, box=SAGUENAY_BOX, beta=BETA_50HZ, z_max_pu=None, y_min_pu=1e-5):
    """Équivalent de Kron du réseau < 315 kV d'une région, ramené aux bus de la région qui existent dans
    le modèle (bus frontières). Les paramètres d'elec_full sont d'abord corrigés : longueurs (Mercator, table T)
    et fréquence (types 50 Hz -> x, b × beta). Retourne la liste des branches équivalentes (pu base 1 MVA)
    et les shunts équivalents. Les charges des bus éliminés sont déjà dans le modèle (réduction de Nathan)."""
    B = full.buses
    inbox = B.index[(B.y.between(*box["lat"])) & (B.x.between(*box["lon"]))]
    model_b = set(model.buses.index)
    region = set(inbox)
    # ordre trié : un `set` de chaînes change d'ordre d'un processus Python à l'autre, ce qui changeait l'ordre
    # des lignes ajoutées et donc la solution (dégénérée) du LOPF -> résultats non reproductibles
    elim = sorted(b for b in region if b not in model_b)
    bound = sorted(b for b in region if b in model_b)
    idx = {b: i for i, b in enumerate(bound + elim)}
    N = len(idx)
    Y = np.zeros((N, N), complex)
    with quiet():
        full.calculate_dependent_values()
    L = full.lines[full.lines.bus0.isin(region) & full.lines.bus1.isin(region)]
    used_lines = []
    for name, r in L.iterrows():
        a, b = r.bus0, r.bus1
        if a in model_b and b in model_b:
            continue  # branche déjà dans le modèle (corrigée par ailleurs)
        f = T.at[name, "factor"] if name in T.index else 1.0
        fb = beta if r.v_nom < 315 else 1.0
        vn2 = r.v_nom ** 2
        z = (r.r * f + 1j * r.x * f * fb) / vn2
        if abs(z) == 0:
            continue
        ysh = 1j * r.b * f * fb * vn2 / 2
        i, j = idx[a], idx[b]
        Y[i, i] += 1 / z + ysh; Y[j, j] += 1 / z + ysh; Y[i, j] -= 1 / z; Y[j, i] -= 1 / z
        used_lines.append(name)
    Tr = full.transformers[full.transformers.bus0.isin(region) & full.transformers.bus1.isin(region)]
    for name, r in Tr.iterrows():
        a, b = r.bus0, r.bus1
        if a in model_b and b in model_b:
            continue
        z = (r.r + 1j * r.x) / r.s_nom
        if abs(z) == 0:
            z = 1j * 0.1 / r.s_nom
        i, j = idx[a], idx[b]
        Y[i, i] += 1 / z; Y[j, j] += 1 / z; Y[i, j] -= 1 / z; Y[j, i] -= 1 / z
    nbd = len(bound)
    Ybb, Ybe, Yeb, Yee = Y[:nbd, :nbd], Y[:nbd, nbd:], Y[nbd:, :nbd], Y[nbd:, nbd:]
    # bus éliminés isolés (aucune branche) : on les retire
    keep = np.where(np.abs(Yee).sum(axis=1) > 0)[0]
    Yee = Yee[np.ix_(keep, keep)]; Ybe = Ybe[:, keep]; Yeb = Yeb[keep, :]
    Yeq = Ybb - Ybe @ np.linalg.solve(Yee, Yeb)
    br = []
    for i in range(nbd):
        for j in range(i + 1, nbd):
            y = -Yeq[i, j]
            if abs(y) > y_min_pu:
                br.append(dict(bus0=bound[i], bus1=bound[j], z_pu=1 / y))
    sh = {bound[i]: Yeq[i].sum() for i in range(nbd)}
    return pd.DataFrame(br), sh, dict(boundary=bound, eliminated=len(elim), lines=len(used_lines))

def add_kron_equivalent(n, br, sh, prefix="KRON", log=None):
    """Ajoute les branches équivalentes : Line entre bus de même tension (x, r en ohm), Transformer sinon
    (x, r en pu sur s_nom = 1 000 MVA). Les shunts équivalents (charge capacitive des lignes éliminées) aussi."""
    added = []
    for k, r in br.iterrows():
        a, b, z = r.bus0, r.bus1, r.z_pu
        va, vb = n.buses.v_nom[a], n.buses.v_nom[b]
        nm = f"{prefix} {a}-{b}"
        with quiet():
            if va == vb:
                n.add("Line", nm, bus0=a, bus1=b, r=max(z.real, 1e-6) * va ** 2, x=z.imag * va ** 2, b=0.0,
                      s_nom=5000.0, length=0.0, num_parallel=1.0, v_nom=va, carrier="AC")
            else:
                hv, lv = (a, b) if va > vb else (b, a)
                n.add("Transformer", nm, bus0=hv, bus1=lv, r=max(z.real, 1e-6) * 1000.0, x=z.imag * 1000.0,
                      s_nom=1000.0, model="t")
        added.append(nm)
    for b, y in sh.items():
        if abs(y.imag) > 1e-6:
            with quiet():
                n.add("ShuntImpedance", f"{prefix} sh {b}", bus=b, b=y.imag / n.buses.v_nom[b] ** 2, g=0.0)
    n.calculate_dependent_values()
    if log is not None:
        log.append(("E2b équivalent de Kron (Saguenay)", f"{len(added)} branches équivalentes"))
    return added



# ============================================================================ E3 circuits (Nathan)
def correct_circuits(n, log=None):
    with quiet():
        diff = _csc.correct_circuits(n)
    if log is not None:
        log.append(("E3 circuits", f"{len(diff)} lignes corrigées (script de Nathan)"))
    return diff



# ============================================================================ E4 compensation série
SERIES_CORRIDORS = {
    "Némiscau-Abitibi (457-312)": ["163", "1655", "2386"],
    "Abitibi-La Vérendrye (312-308)": ["880", "1935", "2410"],
    "Albanel-Chibougamau (345-310, 3319)": ["1653", "2387", "2388"],
    "Chibougamau-Chamouchouane (310-286)": ["1605", "1652", "2272"],
    "Churchill Falls-Montagnais (469-306)": ["1568", "1569", "1570"],
    "Manicouagan-Lévis (114-148, via Bergeronnes)": ["1049", "1448"],
    "Micoua-Laurentides (651-25, via Bergeronnes)": ["3756"],
    "Manicouagan-Québec (114-1291, via Bergeronnes)": ["2408"],
}

SERIES_FRACTION = 0.40

def x_uncompensated(n, line):
    ref = load_hq_params()
    tier = VOLTAGE_TIER_MAP[n.lines.at[line, "v_nom"]]
    return ref.loc[tier, "xl_ohm_per_km"] * n.lines.at[line, "length"] / n.lines.at[line, "num_parallel"]

def series_compensation(n, corridors=None, fraction=SERIES_FRACTION, log=None):
    corridors = SERIES_CORRIDORS if corridors is None else corridors
    done = []
    for corr, lines in corridors.items():
        for l in lines:
            if l in n.lines.index:
                n.lines.at[l, "x"] = x_uncompensated(n, l) * (1 - fraction)
                done.append((corr, l))
    n.calculate_dependent_values()
    if log is not None:
        log.append(("E4 compensation série", f"{fraction:.0%} sur {len(done)} lignes ({len(corridors)} corridors)"))
    return pd.DataFrame(done, columns=["corridor", "line"])



# ============================================================================ E5 capacités + LOPF
ST_CLAIR_MARGIN = 3.0   # même marge que run_lopf_main_island.py

def st_clair_ratings(n, delta_max_deg=None):
    """s_nom = 3 × St Clair (comme Nathan), avec les vraies longueurs. Si delta_max_deg est donné, on
    plafonne aussi par la limite de stabilité statique de chaque ligne : P <= (V^2 / x) sin(delta_max)
    (le LOPF de Nathan ne la voyait pas : la ligne 3360 de la Romaine transitait jusqu'à 90 % de V^2/x)."""
    ref = load_line_params()
    n.lines["s_nom"] = [ST_CLAIR_MARGIN * st_clair_s_nom_mva(r.length, r.v_nom, ref, num_parallel=max(r.num_parallel, 1.0))
                        if r.length > 0 else 5000.0 for _, r in n.lines.iterrows()]
    if delta_max_deg:
        p_stab = n.lines.v_nom ** 2 / n.lines.x.abs() * np.sin(np.radians(delta_max_deg))
        n.lines["s_nom_stclair"] = n.lines["s_nom"]
        n.lines["s_nom"] = np.minimum(n.lines["s_nom"], p_stab)
    n.lines["s_nom_opt"] = n.lines["s_nom"]

def lopf(n, solver="highs", log=None):
    """Relance le LOPF de Nathan (coûts, p_max_pu, délestage déjà dans elec_solved.nc).
    Piège PyPSA 1.3 : un p_set non-NaN fige le dispatch -> on le vide avant optimize()."""
    for c in [n.generators, n.storage_units]:
        c["p_set"] = np.nan
    for t in [n.generators_t, n.storage_units_t]:
        if "p_set" in t and not t["p_set"].empty:
            t["p_set"] = t["p_set"].iloc[:, :0]
    t0 = time.time()
    with quiet():
        status, cond = n.optimize(solver_name=solver, solver_options={"output_flag": False, "log_to_console": False})
    shed = n.generators.index[n.generators.carrier == "load_shedding"]
    shed_mwh = float(n.generators_t.p[shed].sum().sum())
    dem = float(n.loads_t.p_set.sum().sum())
    out = dict(status=status, condition=cond, shed_MWh=shed_mwh, shed_pct=100 * shed_mwh / dem, s=time.time() - t0)
    if log is not None:
        log.append(("E5 LOPF", f"{status}/{cond}, délestage {out['shed_pct']:.3f} %"))
    return out



# ============================================================================ E6 préparation AC
REACTIVE_DEVICES = [
    ("Albanel", "345", "SVC", 2, -105, 300, "Régie R-3712-2009, tableau 2"),
    ("Chamouchouane", "286", "SVC", 2, -300, 300, "Régie R-3712-2009 ; Hitachi Energy (±330)"),
    ("Chibougamau", "310", "SVC", 2, -105, 300, "Régie R-3712-2009"),
    ("Laurentides", "25", "SVC", 1, -100, 350, "Régie R-3712-2009"),
    ("Némiscau", "457", "SVC", 2, -105, 300, "Régie R-3712-2009"),
    ("La Vérendrye", "308", "SVC", 2, -105, 300, "Régie R-3712-2009 ; R-4125-2020 (-100/+300)"),
    ("Lévis (déglaceur)", "148", "SVC", 1, -125, 250, "Régie R-3712-2009"),
    ("Chénier", "136", "SVC", 2, -300, 300, "Hitachi Energy (2011, ±300 chacun)"),
    ("Lévis", "148", "CS", 2, -200, 300, "Régie D-2005-45 (-200/+300 chacun)"),
    ("Duvernay", "195", "CS", 3, -200, 300, "3 CS (Régie R-3876-2014) ; puissance supposée = Lévis"),
    ("Abitibi", "311", "CS", 2, -200, 300, "2 CS (Régie R-3684-2009) ; puissance supposée = Lévis"),
    ("Manicouagan", "114", "CS", 1, -200, 300, "CS (Régie R-3810-2012) ; nombre et puissance supposés"),
]

Q_SYNC = (-0.30, 0.48)   # centrales synchrones (hydro, fil de l'eau, TAG) : Q/p_nom (cos phi 0,9 en surexcitation)

Q_RENEW = (-0.33, 0.33)  # éolien, solaire : cos phi 0,95 en absorption et en production

SYNC = {"hydro", "ror", "OCGT"}

RENEW = {"onwind", "solar", "offwind"}

def reactive_devices_table():
    return pd.DataFrame(REACTIVE_DEVICES, columns=["poste", "bus", "type", "nombre", "qmin_unit", "qmax_unit", "source"])

def bus_q_capability(n, devices=True, sync=Q_SYNC, renew=Q_RENEW):
    """Capacité réactive par bus : centrales (p_nom installée) + stockage hydro + SVC/CS réels."""
    g = n.generators[(n.generators.carrier != "load_shedding") & ~n.generators.index.str.contains("placeholder")]
    rows = []
    for _, r in g.iterrows():
        if r.carrier in SYNC:
            rows.append((r.bus, "centrale", sync[0] * r.p_nom, sync[1] * r.p_nom))
        elif r.carrier in RENEW:
            rows.append((r.bus, "éolien/solaire", renew[0] * r.p_nom, renew[1] * r.p_nom))
    for _, r in n.storage_units.iterrows():
        rows.append((r.bus, "centrale", sync[0] * r.p_nom, sync[1] * r.p_nom))
    if devices:
        for _, d in reactive_devices_table().iterrows():
            if d.bus in n.buses.index:
                rows.append((d.bus, d.type, d.nombre * d.qmin_unit, d.nombre * d.qmax_unit))
    C = pd.DataFrame(rows, columns=["bus", "source", "qmin", "qmax"])
    return C

def add_line_reactors(n, min_vnom=700.0, min_km_315=100.0):
    """Une inductance shunt par bus, égale à la moitié de la susceptance de chaque ligne qui y arrive
    (k = 1) : toutes les lignes 735/765 kV, et les lignes 315 kV de plus de min_km_315 km (lignes longues
    peu chargées, ex. Fermont). Le taux k réellement enclenché est choisi heure par heure (b = k × b_unit)."""
    L = n.lines[(n.lines.v_nom >= min_vnom) | (n.lines.length >= min_km_315)]
    b_end = pd.concat([L.set_index("bus0").b, L.set_index("bus1").b]).groupby(level=0).sum() / 2.0
    names = [f"REACT {b}" for b in b_end.index]
    with quiet():
        n.add("ShuntImpedance", names, bus=b_end.index.values, b=0.0, g=0.0)
    n.shunt_impedances.loc[names, "b_unit"] = -b_end.values
    n.shunt_impedances.loc[names, "kind"] = "reactor"
    return float((b_end * (n.buses.v_nom[b_end.index] * 1e3) ** 2).sum() / 1e6)  # MVAr à k = 1, 1 pu

def add_capacitor_banks(n, banks):
    """Bancs de condensateurs commutables aux postes de charge : banks = {bus: MVAr à 1 pu}."""
    for b, q in banks.items():
        nm = f"CAP {b}"
        bu = q * 1e6 / (n.buses.v_nom[b] * 1e3) ** 2
        if nm in n.shunt_impedances.index:
            n.shunt_impedances.at[nm, "b_unit"] += bu
        else:
            with quiet():
                n.add("ShuntImpedance", nm, bus=b, b=0.0, g=0.0)
            n.shunt_impedances.loc[nm, "b_unit"] = bu
            n.shunt_impedances.loc[nm, "kind"] = "capacitor"

def capacitor_banks_from_load(n, alpha=0.6, min_q=20.0, step=None):
    """Bancs de condensateurs aux postes de charge : alpha × Q de la charge locale à la pointe.
    Représente la compensation de la distribution et du réseau régional (non modélisée), qui ramène
    le facteur de puissance vu du 315 kV vers 0,98-1,0. alpha = 0,57 -> cos phi 0,99 à la pointe."""
    lq = n.loads_t.q_set.T.groupby(n.loads.bus).sum().max(axis=1)
    lq = lq[lq.index.isin(n.buses.index[n.buses.carrier == "AC"])]
    banks = (alpha * lq)[alpha * lq >= min_q]
    if step:
        banks = (np.ceil(banks / step) * step)
    return banks.to_dict()

def prepare_ac(n, load_pf=0.95, v_set=1.0, devices=True, log=None):
    """Prépare le réseau (issu du LOPF) pour l'AC :
    1. préparation de Nathan (dispatch -> p_set, r des transfos, Q des charges) ;
    2. Q des charges recalculé avec load_pf ;
    3. un régulateur par bus avec capacité Q : générateur fictif « REG bus » (P = 0, PV, limites Q =
       somme des capacités du bus) ; les vraies centrales passent en PQ à Q = 0 (elles gardent leur P) ;
       le slack de Nathan reste le slack (son Q est vérifié a posteriori) ;
    4. inductances de ligne 735 kV (commutables, k appliqué heure par heure)."""
    with quiet():
        run_pf.use_lopf_dispatch_as_setpoints(n)
        n.determine_network_topology()
        run_pf.prepare_for_ac_pf(n, pv_min_capacity=0.0, pv_max_load_ratio=1.0)
    tan = np.tan(np.arccos(load_pf))
    n.loads_t.q_set = n.loads_t.p_set * tan
    n.loads["q_set"] = n.loads["p_set"] * tan
    # placeholders de Nathan inutiles ici (le régulateur prend leur rôle)
    ph = n.generators.index[n.generators.index.str.contains("storage-pv-placeholder")]
    with quiet():
        n.remove("Generator", ph)
    slack = list(n.generators.index[n.generators.control == "Slack"])
    real = n.generators.index[n.generators.carrier != "load_shedding"].difference(slack)
    n.generators.loc[real, "control"] = "PQ"
    n.generators.loc[real, "q_set"] = 0.0
    if not n.generators_t.q_set.empty:
        n.generators_t.q_set = n.generators_t.q_set.iloc[:, :0]
    C = bus_q_capability(n, devices=devices)
    cap = C.groupby("bus")[["qmin", "qmax"]].sum()
    slack_buses = set(n.generators.bus[slack])
    regs = []
    for b, r in cap.iterrows():
        if b in slack_buses:
            continue
        nm = f"REG {b}"
        with quiet():
            n.add("Generator", nm, bus=b, carrier="regulator", p_nom=0.0, p_set=0.0, control="PV")
        regs.append(nm)
    n.generators.loc[regs, "q_max"] = cap.loc[[r[4:] for r in regs], "qmax"].values
    n.generators.loc[regs, "q_min"] = cap.loc[[r[4:] for r in regs], "qmin"].values
    for s in slack:
        b = n.generators.at[s, "bus"]
        n.generators.loc[s, ["q_min", "q_max"]] = cap.loc[b, ["qmin", "qmax"]].values if b in cap.index else (np.nan, np.nan)
    n.generators_t.q_set = pd.DataFrame(0.0, index=n.snapshots, columns=regs)
    n.buses["v_mag_pu_set"] = 1.0
    n.buses.loc[cap.index, "v_mag_pu_set"] = v_set
    q_react = add_line_reactors(n)
    if log is not None:
        log.append(("E6 préparation AC", f"{len(regs)} régulateurs, cos phi charges {load_pf}, "
                    f"inductances de ligne 735 kV : {q_react:.0f} MVAr à k = 1"))
    return dict(regulators=regs, slack=slack, capability=C, reactors_MVAr_k1=q_react)



# ============================================================================ pipeline E1 -> E6
def build(source=None, full=None, T=None, series_fraction=SERIES_FRACTION, corridors=None, circuits=True,
          churchill=True, kron=True, lengths=True, rerun_lopf=True, load_pf=0.95, v_set=1.02, devices=True,
          alpha_banks=1.0, delta_max_deg=40.0, log=None):
    """E1 -> E6 avec les choix du cas final par défaut (chaque argument permet de retirer ou changer une
    hypothèse, cf. sensibilites.py). Retourne (n, info) ; info["log"] résume chaque étape.
    T = table Mercator (calculée depuis elec_full si absente)."""
    log = [] if log is None else log
    n = load(source or os.path.join(CURRENT, "elec_solved.nc"))
    info = {}
    if lengths:
        if T is None:
            T = mercator_table(load(full or os.path.join(CURRENT, "elec_full.nc")))
        info["lengths"] = fix_lengths(n, T, log)
    if churchill:
        fix_churchill_falls(n, log)
    if kron:
        full_n = load(full or os.path.join(CURRENT, "elec_full.nc"))
        if T is None:
            T = mercator_table(full_n)
        regions = KRON_REGIONS if kron is True else {k: KRON_REGIONS[k] for k in kron}
        info["kron"] = {}
        for reg, box in regions.items():
            br, sh, kinfo = kron_equivalent(full_n, n, T, box=box)
            info["kron"][reg] = (br, sh, kinfo)
            add_kron_equivalent(n, br, sh, prefix=f"KRON {reg}")
            log.append((f"E2b équivalent de Kron ({reg})", f"{kinfo['eliminated']} bus < 315 kV éliminés -> {len(br)} branches équivalentes"))
    if circuits:
        info["circuits"] = correct_circuits(n, log)
    if series_fraction:
        info["series"] = series_compensation(n, corridors, series_fraction, log)
    if rerun_lopf:
        st_clair_ratings(n, delta_max_deg)
        info["lopf"] = lopf(n, log=log)
    info["ac"] = prepare_ac(n, load_pf=load_pf, v_set=v_set, devices=devices, log=log)
    if alpha_banks:
        banks = capacitor_banks_from_load(n, alpha=alpha_banks)
        add_capacitor_banks(n, banks)
        info["banks"] = banks
        log.append(("E6c condensateurs", f"{len(banks)} postes de charge, {sum(banks.values()):.0f} MVAr installés (alpha = {alpha_banks})"))
    info["log"] = log
    return n, info



# ============================================================================ E7 moteur Newton-Raphson et exploitation heure par heure
K_PEAK, K_LIGHT = 0.37, 0.85   # HQ 1996 : 11 000 / 30 000 Mvar enclenchés à la pointe, 25 500 / 30 000 au total

class ACEngine:
    """Newton-Raphson polaire (comme nr.py du cours), sur la matrice Y que PyPSA construit pour le
    sous-réseau AC principal (même modèle, même convention : MW / MVAr, Y en pu sur base 1 MVA).
    Les shunts commutables (inductances, condensateurs) sont ajoutés sur la diagonale heure par heure.
    But : quelques ms par résolution, pour pouvoir commuter limites Q et inductances heure par heure.
    Les résultats finaux sont ensuite revérifiés par le solveur de PyPSA (validate_with_pypsa)."""

    def __init__(self, n, regs):
        with quiet():
            n.determine_network_topology()
        self.n = n
        subs = list(n.sub_networks.obj)
        self.sub = max(subs, key=lambda s: len(s.buses_i()))
        sh = n.shunt_impedances
        self.sw = sh.index[sh.kind.isin(["reactor", "capacitor"])] if "kind" in sh else sh.index[:0]
        b_save = sh.loc[self.sw, "b"].copy()
        sh.loc[self.sw, "b"] = 0.0
        with quiet():
            n.calculate_dependent_values()
            self.sub.find_bus_controls()
            self.sub.calculate_Y()
        sh.loc[self.sw, "b"] = b_save
        self.buses = list(self.sub.buses_o)
        self.ix = {b: i for i, b in enumerate(self.buses)}
        self.nb = len(self.buses)
        self.Y0 = self.sub.Y.toarray()
        self.slack_bus = self.sub.slack_bus
        self.sl = self.ix[self.slack_bus]
        vn = n.buses.v_nom.reindex(self.buses).values
        self.vnom = vn
        # susceptances unitaires des shunts commutables, en pu base 1 MVA : b_pu = b [S] * v_nom[kV]^2
        self.b_react = np.zeros(self.nb); self.b_cap = np.zeros(self.nb)
        for s in self.sw:
            b = sh.at[s, "bus"]
            if b in self.ix:
                bp = sh.at[s, "b_unit"] * n.buses.v_nom[b] ** 2
                (self.b_react if sh.at[s, "kind"] == "reactor" else self.b_cap)[self.ix[b]] += bp
        self.regs = [g for g in regs if n.generators.at[g, "bus"] in self.ix]
        self.reg_bus = np.array([self.ix[n.generators.at[g, "bus"]] for g in self.regs])
        self.qmin = n.generators.loc[self.regs, "q_min"].values.astype(float)
        self.qmax = n.generators.loc[self.regs, "q_max"].values.astype(float)
        self.vset = n.buses.v_mag_pu_set.reindex(self.buses).values.astype(float)
        self._prep_injections()

    def _prep_injections(self):
        n, ix, sn = self.n, self.ix, self.n.snapshots
        P = pd.DataFrame(0.0, index=sn, columns=self.buses)
        Q = pd.DataFrame(0.0, index=sn, columns=self.buses)

        def add(df, ts, buses, sign=1.0):
            for c in ts.columns:
                b = buses[c]
                if b in ix:
                    df[b] += sign * ts[c].values

        g = n.generators
        gp = n.generators_t.p_set.reindex(columns=g.index)
        gp = gp.fillna(pd.Series(g.p_set.fillna(0.0).values, index=g.index)).fillna(0.0)
        reg_or_slack = set(self.regs) | set(g.index[g.control == "Slack"])
        add(P, gp[[c for c in gp.columns if c not in reg_or_slack]], g.bus)
        # Q fixe des générateurs non régulateurs (q_set statique ; 0 pour les vraies centrales ici)
        gq = pd.DataFrame({c: np.full(len(sn), float(g.at[c, "q_set"] or 0.0)) for c in g.index if c not in reg_or_slack}, index=sn)
        add(Q, gq, g.bus)
        su = n.storage_units
        if len(su):
            sp = n.storage_units_t.p_set.reindex(columns=su.index).fillna(0.0)
            add(P, sp, su.bus)
        lp = n.loads_t.p_set.reindex(columns=n.loads.index).fillna(pd.Series(n.loads.p_set.values, index=n.loads.index))
        lq = n.loads_t.q_set.reindex(columns=n.loads.index).fillna(pd.Series(n.loads.q_set.values, index=n.loads.index))
        add(P, lp, n.loads.bus, -1.0); add(Q, lq, n.loads.bus, -1.0)
        if len(n.links):
            p0 = n.links_t.p_set.reindex(columns=n.links.index).fillna(0.0)
            add(P, p0, n.links.bus0, -1.0)
            add(P, p0 * n.links.efficiency.values, n.links.bus1, +1.0)
        self.P, self.Q = P, Q

    # ---------------------------------------------------------------- Newton
    def _dSdV(self, Y, V):
        I = Y @ V
        dV = np.diag(V); dI = np.diag(I); dVn = np.diag(V / np.abs(V))
        dth = 1j * dV @ np.conj(dI - Y @ dV)
        dvm = dV @ np.conj(Y @ dVn) + np.conj(dI) @ dVn
        return dth, dvm

    def dc_start(self, Y, P):
        B = -Y.imag.copy()
        np.fill_diagonal(B, 0.0); np.fill_diagonal(B, -B.sum(axis=1))
        ns = [i for i in range(self.nb) if i != self.sl]
        th = np.zeros(self.nb)
        th[ns] = np.linalg.solve(-B[np.ix_(ns, ns)], P[ns])
        return th

    def newton(self, Y, S, V0, pv, pq, tol=1e-3, maxit=30):
        V = V0.copy()
        pvpq = np.r_[pv, pq].astype(int)
        npv = len(pvpq)
        F = None
        for it in range(maxit + 1):
            mis = V * np.conj(Y @ V) - S
            F = np.r_[mis[pvpq].real, mis[pq].imag]
            err = np.abs(F).max()
            if not np.isfinite(err):
                return V, False, it
            if err < tol:
                return V, True, it
            dth, dvm = self._dSdV(Y, V)
            J = np.block([[dth[np.ix_(pvpq, pvpq)].real, dvm[np.ix_(pvpq, pq)].real],
                          [dth[np.ix_(pq, pvpq)].imag, dvm[np.ix_(pq, pq)].imag]])
            try:
                dx = np.linalg.solve(J, -F)
            except np.linalg.LinAlgError:
                return V, False, it
            # pas amorti (backtracking) si la norme du résidu augmente
            th0, vm0 = np.angle(V), np.abs(V)
            for a in (1.0, 0.5, 0.25, 0.1):
                th = th0.copy(); vm = vm0.copy()
                th[pvpq] += a * dx[:npv]; vm[pq] += a * dx[npv:]
                Vt = vm * np.exp(1j * th)
                mt = Vt * np.conj(Y @ Vt) - S
                Ft = np.r_[mt[pvpq].real, mt[pq].imag]
                if np.isfinite(Ft).all() and np.abs(Ft).max() < err or a == 0.1:
                    V = Vt; break
            if np.abs(V).min() < 0.3 or np.abs(V).max() > 2.0:
                return V, False, it
        return V, False, maxit

    # ---------------------------------------------------------------- une heure
    def solve(self, sn, k=0.5, c=0.0, max_rounds=40, Vstart=None, back_switch=True, max_flips=2):
        """AC PF d'une heure avec limites Q (bascule PV -> PQ, et retour PQ -> PV si la tension le permet).
        Un régulateur qui a déjà changé d'état max_flips fois est bloqué (évite les oscillations).
        Retourne dict(V, ok, q_reg, state, n_viol, ...) ; `state` est l'état EXACT utilisé pour calculer V,
        et n_viol le nombre de régulateurs PV encore hors limites dans cette solution (0 si tout est respecté)."""
        kv = np.broadcast_to(np.asarray(k, float), (self.nb,))
        cv = np.broadcast_to(np.asarray(c, float), (self.nb,))
        Y = self.Y0 + np.diag(1j * (kv * self.b_react + cv * self.b_cap))
        P = self.P.loc[sn].values; Qf = self.Q.loc[sn].values
        nreg = len(self.regs)
        state = np.zeros(nreg, int)   # 0 = PV, +1 = PQ à Qmax, -1 = PQ à Qmin
        flips = np.zeros(nreg, int)
        if Vstart is None:
            th = self.dc_start(Y, P)
            vm = np.ones(self.nb)
            vm[self.reg_bus] = self.vset[self.reg_bus]; vm[self.sl] = self.vset[self.sl]
            V = vm * np.exp(1j * th)
        else:
            V = Vstart.copy()
        ok = False; nsw = 0; state_v = state.copy(); V_ok = None
        for rnd in range(max_rounds):
            pv_set = {self.reg_bus[j] for j in range(nreg) if state[j] == 0}
            pv_b = sorted(pv_set)
            pq_b = [i for i in range(self.nb) if i != self.sl and i not in pv_set]
            Qs = Qf.copy()
            for j in np.where(state != 0)[0]:
                Qs[self.reg_bus[j]] += self.qmax[j] if state[j] > 0 else self.qmin[j]
            Vm = np.abs(V); Vm[pv_b] = self.vset[pv_b]; Vm[self.sl] = self.vset[self.sl]
            V = Vm * np.exp(1j * np.angle(V))
            V, ok, it = self.newton(Y, P + 1j * Qs, V, np.array(pv_b, int), np.array(pq_b, int))
            if not ok:
                break
            V_ok, state_v = V.copy(), state.copy()     # solution et état cohérents
            qbus = (V * np.conj(Y @ V)).imag - Qf      # Q fourni par le régulateur du bus (MVAr)
            changed = False
            viol = []
            for j in range(nreg):
                b = self.reg_bus[j]
                if state[j] == 0:
                    if qbus[b] > self.qmax[j] + 1.0:
                        viol.append((qbus[b] - self.qmax[j], j, +1))
                    elif qbus[b] < self.qmin[j] - 1.0:
                        viol.append((self.qmin[j] - qbus[b], j, -1))
                elif back_switch and flips[j] < max_flips:
                    vmb = abs(V[b])
                    if (state[j] > 0 and vmb > self.vset[b] + 1e-4) or (state[j] < 0 and vmb < self.vset[b] - 1e-4):
                        state[j] = 0; flips[j] += 1; changed = True
            for _, j, sgn in sorted(viol, reverse=True)[:3]:
                state[j] = sgn; flips[j] += 1; changed = True; nsw += 1
            if not changed:
                break
        if V_ok is None:
            return dict(V=V, ok=False, q_reg=np.full(nreg, np.nan), state=state_v, n_switch=nsw, n_viol=-1,
                        p_slack=np.nan, q_slack=np.nan, q_react=np.zeros(self.nb), q_cap=np.zeros(self.nb))
        V = V_ok
        Sc = V * np.conj(Y @ V)
        qbus = Sc.imag - Qf
        q_reg = np.array([qbus[b] if state_v[j] == 0 else (self.qmax[j] if state_v[j] > 0 else self.qmin[j])
                          for j, b in enumerate(self.reg_bus)])
        n_viol = int(sum(1 for j in range(nreg) if state_v[j] == 0 and (q_reg[j] > self.qmax[j] + 1.0 or q_reg[j] < self.qmin[j] - 1.0)))
        return dict(V=V, ok=True, q_reg=q_reg, state=state_v, n_switch=nsw, n_viol=n_viol,
                    p_slack=Sc[self.sl].real - P[self.sl], q_slack=Sc[self.sl].imag - Qf[self.sl],
                    q_react=-kv * self.b_react * np.abs(V) ** 2, q_cap=cv * self.b_cap * np.abs(V) ** 2)

def neighbours(eng):
    """Bus voisins (lignes et transformateurs) dans l'ordre du moteur."""
    n = eng.n; nb = {i: set() for i in range(eng.nb)}
    for df in (n.lines, n.transformers):
        for a, b in zip(df.bus0, df.bus1):
            if a in eng.ix and b in eng.ix:
                nb[eng.ix[a]].add(eng.ix[b]); nb[eng.ix[b]].add(eng.ix[a])
    return nb

def solve_hour_reserve(eng, sn, k0, c0, kmax=0.85, kmin=0.0, dk=0.05, dc=0.25, q_band=0.7,
                       v_lo=0.98, v_hi=1.04, cap_lo=0.97, cap_hi=1.03, iters_a=25, iters_b=12, nbrs=None):
    """Une heure, en deux phases (logique d'exploitation d'HQ : les inductances et condensateurs sont
    manoeuvrés pour garder de la réserve dans les SVC / compensateurs synchrones et les centrales) :
    A. limites Q relâchées (Newton converge toujours mieux) : chaque poste qui a des shunts lit soit le Q
       de son régulateur (s'il en a un, ou s'il en a un voisin), soit sa tension ; si le régulateur produit
       plus de q_band × Qmax -> on retire des inductances / on enclenche des condensateurs ; s'il absorbe plus
       de q_band × |Qmin| -> l'inverse ; sans régulateur, même logique sur la tension (v_lo, v_hi).
    B. vraies limites Q (bascule PV -> PQ), puis ajustement local sur la tension (cap_lo/cap_hi, v_lo/v_hi)."""
    nbrs = nbrs or neighbours(eng)
    hasR = eng.b_react != 0; hasC = eng.b_cap != 0
    k = np.where(hasR, k0, 0.0).astype(float); c = np.where(hasC, c0, 0.0).astype(float)
    reg_of = {}
    for j, b in enumerate(eng.reg_bus):
        reg_of.setdefault(b, []).append(j)
    qmin0, qmax0 = eng.qmin.copy(), eng.qmax.copy()
    # ---------------- phase A : limites relâchées
    eng.qmin[:] = -1e12; eng.qmax[:] = 1e12
    V = None
    try:
        for ia in range(iters_a):
            r = eng.solve(sn, k=k, c=c, Vstart=V)
            if not r["ok"]:
                V = None
                k[hasR] = np.maximum(kmin, k[hasR] - dk); c[hasC] = np.minimum(1.0, c[hasC] + dc)
                continue
            V = r["V"]; vm = np.abs(V); q = r["q_reg"]
            # demande de réactif de chaque régulateur : +1 veut plus de Q, -1 en absorbe trop, 0 ok
            dem = np.zeros(eng.nb)
            for j, b in enumerate(eng.reg_bus):
                if q[j] > q_band * qmax0[j] and qmax0[j] > 0:
                    dem[b] += (q[j] - q_band * qmax0[j])
                elif q[j] < q_band * qmin0[j] and qmin0[j] < 0:
                    dem[b] -= (q_band * qmin0[j] - q[j])
            kn, cn = k.copy(), c.copy()
            for s in np.where(hasR | hasC)[0]:
                if s in reg_of:
                    sig = dem[s]
                else:
                    sig = sum(dem[t] for t in nbrs[s] if t in reg_of)
                    if sig == 0:
                        sig = 1.0 if vm[s] < v_lo else (-1.0 if vm[s] > v_hi else 0.0)
                if sig > 0:      # il manque du Q ici
                    if hasR[s] and k[s] > kmin + 1e-9: kn[s] = max(kmin, k[s] - dk)
                    elif hasC[s] and c[s] < 1 - 1e-9: cn[s] = min(1.0, c[s] + dc)
                elif sig < 0:    # trop de Q ici
                    if hasC[s] and c[s] > 1e-9: cn[s] = max(0.0, c[s] - dc)
                    elif hasR[s] and k[s] < kmax - 1e-9: kn[s] = min(kmax, k[s] + dk)
            if np.allclose(kn, k) and np.allclose(cn, c):
                break
            k, c = kn, cn
    finally:
        eng.qmin[:], eng.qmax[:] = qmin0, qmax0
    # ---------------- phase B : vraies limites
    last = None
    for ib in range(iters_b):
        r = eng.solve(sn, k=k, c=c)
        if not r["ok"]:
            if last is not None:
                break
            k[hasR] = np.maximum(kmin, k[hasR] - dk); c[hasC] = np.minimum(1.0, c[hasC] + dc)
            continue
        if last is None or r["n_viol"] <= last["n_viol"]:
            last = dict(r, k=k.copy(), c=c.copy())
        vm = np.abs(r["V"])
        kn, cn = k.copy(), c.copy()
        up = hasR & (vm > v_hi); dn = hasR & (vm < v_lo)
        kn[up] = np.minimum(kmax, kn[up] + dk); kn[dn] = np.maximum(kmin, kn[dn] - dk)
        cu = hasC & (vm < cap_lo); cd = hasC & (vm > cap_hi)
        cn[cu] = np.minimum(1.0, cn[cu] + dc); cn[cd] = np.maximum(0.0, cn[cd] - dc)
        if np.allclose(kn, k) and np.allclose(cn, c):
            break
        k, c = kn, cn
    if last is None:
        return dict(ok=False, k=k, c=c, iters=(ia + 1, ib + 1))
    last["iters"] = (ia + 1, ib + 1)
    return last

def solve_week_reserve(eng, snapshots=None, k_peak=K_PEAK, k_light=K_LIGHT, v_hard=(0.90, 1.10), band=(0.95, 1.05),
                       verbose=False, **kw):
    n = eng.n
    snapshots = list(n.snapshots if snapshots is None else snapshots)
    load = n.loads_t.p_set.sum(axis=1)
    lam = (load - load.min()) / (load.max() - load.min())
    b735 = np.array([i for i, b in enumerate(eng.buses) if n.buses.v_nom[b] >= 700])
    nbrs = neighbours(eng)
    rows, Vs, Qr, Ks, Cs, St = [], {}, {}, {}, {}, {}
    t0 = time.time()
    for sn in snapshots:
        k0 = k_light - (k_light - k_peak) * lam[sn]
        c0 = float(np.clip(lam[sn], 0, 1))
        r = solve_hour_reserve(eng, sn, k0, c0, nbrs=nbrs, **kw)
        if not r["ok"]:
            rows.append(dict(snapshot=sn, converged=False, newton_ok=False, iters=str(r["iters"])))
            if verbose: print(sn, "DIVERGE", f"{time.time() - t0:.0f}s", flush=True)
            continue
        vm = np.abs(r["V"]); q_reg = r["q_reg"]
        sat = int(((q_reg >= eng.qmax - 1) | (q_reg <= eng.qmin + 1)).sum())
        # « convergé » = Newton converge ET tous les régulateurs restés en PV sont dans leurs limites Q
        rows.append(dict(snapshot=sn, converged=bool(r["n_viol"] == 0), newton_ok=True, iters=str(r["iters"]),
                         n_pq=int((r["state"] != 0).sum()), n_sat=sat, regulateurs_hors_limites=int(r["n_viol"]),
                         v_min=vm.min(), v_max=vm.max(), v735_min=vm[b735].min(), v735_max=vm[b735].max(),
                         p_slack=r["p_slack"], q_slack=r["q_slack"],
                         q_react_MVAr=float(-(r["k"] * eng.b_react * vm ** 2).sum()),
                         q_cap_MVAr=float((r["c"] * eng.b_cap * vm ** 2).sum()),
                         k_eq=float((r["k"] * eng.b_react).sum() / eng.b_react.sum()),
                         c_eq=float((r["c"] * eng.b_cap).sum() / eng.b_cap.sum()) if eng.b_cap.sum() else 0.0,
                         in_hard=bool((vm >= v_hard[0]).all() and (vm <= v_hard[1]).all()),
                         in_band=bool((vm >= band[0]).all() and (vm <= band[1]).all())))
        Vs[sn] = r["V"]; Qr[sn] = q_reg; Ks[sn] = r["k"]; Cs[sn] = r["c"]; St[sn] = r["state"]
        if verbose: print(sn, "ok", round(vm.min(), 3), round(vm.max(), 3), f"{time.time() - t0:.0f}s", flush=True)
    T = pd.DataFrame(rows).set_index("snapshot")
    mk = lambda d, cols: (pd.DataFrame(d).T.set_axis(cols, axis=1) if d else pd.DataFrame(columns=cols))
    return T, mk(Vs, eng.buses), mk(Qr, eng.regs), mk(Ks, eng.buses), mk(Cs, eng.buses), mk(St, eng.regs)



# ============================================================================ validation par le solveur de PyPSA
def complementarity_violations(eng, r, tol=0.01):
    """Régulateurs en butée dont la tension est du « mauvais » côté de la consigne (à Qmax avec V > consigne,
    ou à Qmin avec V < consigne) de plus de `tol` pu : physiquement, ils devraient revenir en PV. MATPOWER et
    pandapower ne revérifient pas ce point ; on le compte pour juger une solution près du nez de la courbe P-V."""
    V = np.abs(r["V"]); out = []
    for j, b in enumerate(eng.reg_bus):
        s = r["state"][j]
        if (s > 0 and V[b] > eng.vset[b] + tol) or (s < 0 and V[b] < eng.vset[b] - tol):
            out.append((eng.regs[j], float(V[b] - eng.vset[b])))
    return out


def load_margin(eng, sn, k0=None, c0=None, mode="remanoeuvre", step=0.005, lam_max=1.3, compl_tol=0.01):
    """Marge de charge d'une heure (courbe P-V) : charge ET production multipliées par lambda, par pas
    de `step`, jusqu'au premier pas sans solution qui respecte les limites Q (Newton converge, aucun
    régulateur resté en PV hors de ses limites, et aucun régulateur en butée du mauvais côté de sa consigne). mode = "fige" (shunts de l'heure gardés, départ DC à
    chaque pas) ou "remanoeuvre" (shunts re-commandés à chaque pas par solve_hour_reserve).
    Retourne un DataFrame (lam, v_min, bus_v_min, n_pq) des pas réussis ; eng.P/eng.Q sont restaurés."""
    P0, Q0 = eng.P.loc[sn].copy(), eng.Q.loc[sn].copy()
    rows = []
    try:
        for lam in np.round(np.arange(1.0, lam_max + 1e-9, step), 4):
            eng.P.loc[sn] = P0 * lam; eng.Q.loc[sn] = Q0 * lam
            if mode == "fige":
                r = eng.solve(sn, k=k0, c=c0)
            else:
                r = solve_hour_reserve(eng, sn, max(0.0, K_PEAK - (lam - 1)), 1.0)
            if not (r["ok"] and r.get("n_viol", 0) == 0) or complementarity_violations(eng, r, compl_tol):
                break
            vm = np.abs(r["V"])
            rows.append(dict(lam=lam, v_min=vm.min(), bus_v_min=eng.buses[vm.argmin()], n_pq=int((r["state"] != 0).sum())))
    finally:
        eng.P.loc[sn], eng.Q.loc[sn] = P0, Q0
    return pd.DataFrame(rows)


def apply_hour_state(n, eng, sn, K, C, S):
    """Recopie dans le réseau PyPSA l'état final d'une heure : b des shunts commutables, régulateurs
    en PQ (Q à la limite) ou PV."""
    sh = n.shunt_impedances
    kk, cc = K.loc[sn], C.loc[sn]
    for s in eng.sw:
        b = sh.at[s, "bus"]
        f = kk.get(b, 0.0) if sh.at[s, "kind"] == "reactor" else cc.get(b, 0.0)
        sh.at[s, "b"] = f * sh.at[s, "b_unit"]
    st = S.loc[sn]
    for j, g in enumerate(eng.regs):
        s = int(st[g])
        if s == 0:
            n.generators.at[g, "control"] = "PV"
        else:
            n.generators.at[g, "control"] = "PQ"
            n.generators_t.q_set.loc[sn, g] = eng.qmax[j] if s > 0 else eng.qmin[j]

def validate_with_pypsa(n, eng, W, snapshots=None, verbose=False):
    """Pour chaque heure : état final -> n.pf() de PyPSA depuis un départ DC (à froid, sans la solution
    du moteur) -> écart de tension avec le moteur. Écrit aussi les résultats complets dans n (flux, Q...)."""
    T, V, K, C, S = W["T"], W["V"], W["K"], W["C"], W["S"]
    snapshots = [s for s in (snapshots or list(T.index[T.converged.astype(bool)]))]
    rows = []
    t0 = time.time()
    for sn in snapshots:
        apply_hour_state(n, eng, sn, K, C, S)
        with quiet():
            n.lpf([sn])
            r = n.pf([sn], use_seed=True, x_tol=1e-6)
        ok = bool(r.converged.all(axis=1).iloc[0])
        vps = n.buses_t.v_mag_pu.loc[sn, eng.buses].values
        aps = n.buses_t.v_ang.loc[sn, eng.buses].values
        ve = V.loc[sn, eng.buses].values.astype(complex)
        row = dict(snapshot=sn, pypsa_converge=ok, ecart_V_max=float(np.abs(vps - np.abs(ve)).max()),
                   ecart_angle_deg_max=float(np.degrees(np.abs(np.angle(np.exp(1j * (aps - np.angle(ve)))))).max()))
        if row["ecart_V_max"] > 1e-4:
            # près du nez de la courbe P-V, il peut exister deux solutions proches : départ à chaud depuis
            # la solution du moteur -> si PyPSA y reste, c'est bien une solution de ses propres équations
            n.buses_t.v_mag_pu.loc[sn, eng.buses] = np.abs(ve)
            n.buses_t.v_ang.loc[sn, eng.buses] = np.angle(ve)
            with quiet():
                r2 = n.pf([sn], use_seed=True, x_tol=1e-6)
            vps = n.buses_t.v_mag_pu.loc[sn, eng.buses].values
            row.update(depart_chaud_converge=bool(r2.converged.all(axis=1).iloc[0]),
                       ecart_V_max_depart_chaud=float(np.abs(vps - np.abs(ve)).max()))
        rows.append(row)
        if verbose and len(rows) % 24 == 0:
            print(len(rows), f"{time.time() - t0:.0f}s", flush=True)
    return pd.DataFrame(rows).set_index("snapshot")



# ============================================================================ export MATPOWER (AC-CFM) et pandapower
GENTYPE = {"ror": "HY", "hydro": "HY", "AC": "HY", "OCGT": "GT", "onwind": "WN", "solar": "PV", "regulator": "SC"}


def to_ppc(n, eng, W, sn, base_mva=100.0):
    """Cas MATPOWER (dict ppc) d'une heure du sous-réseau AC principal, dans l'état final de l'heure :
    shunts commutés (Bs), régulateurs (générateurs P = 0 avec leurs limites Q et Vg), centrales réelles
    (P du LOPF, Q = 0), slack. Les liaisons HVDC vers les îlots sont des injections fixes (Pd).
    Champ supplémentaire `gentype` (HY, GT, WN, PV, SC) pour AC-CFM ; `bus_names` = bus PyPSA."""
    buses = eng.buses; bi = {b: i + 1 for i, b in enumerate(buses)}
    vn = n.buses.v_nom.reindex(buses)
    T, K, C, S = W["T"], W["K"], W["C"], W["S"]
    # injections fixes hors générateurs (charges, liaisons) -> Pd, Qd
    Pfix = eng.P.loc[sn].copy(); Qfix = eng.Q.loc[sn].copy()
    g = n.generators
    gp = n.generators_t.p_set.reindex(columns=g.index).loc[sn].fillna(pd.Series(g.p_set.fillna(0.0)))
    gens_rows, gtype = [], []
    for name, r in g.iterrows():
        if r.bus not in bi or r.carrier == "load_shedding":
            continue
        p = float(gp.get(name, 0.0) or 0.0)
        if name in eng.regs or r.control == "Slack":
            j = eng.regs.index(name) if name in eng.regs else None
            qmax = eng.qmax[j] if j is not None else (r.q_max if np.isfinite(r.get("q_max", np.nan)) else 9999.0)
            qmin = eng.qmin[j] if j is not None else (r.q_min if np.isfinite(r.get("q_min", np.nan)) else -9999.0)
            vg = float(n.buses.v_mag_pu_set[r.bus])
            gens_rows.append([bi[r.bus], p, 0.0, qmax, qmin, vg, base_mva, 1, max(r.p_nom, p, 1.0), 0.0] + [0] * 11)
            gtype.append("SC" if name in eng.regs else "HY")
        else:
            Pfix[r.bus] -= p   # on la sort de l'injection fixe pour la remettre comme générateur
            gens_rows.append([bi[r.bus], p, 0.0, 0.0, 0.0, float(n.buses.v_mag_pu_set[r.bus]), base_mva, 1, max(r.p_nom, p), 0.0] + [0] * 11)
            gtype.append(GENTYPE.get(r.carrier, "HY"))
    for name, r in n.storage_units.iterrows():
        if r.bus not in bi:
            continue
        p = float(n.storage_units_t.p_set.loc[sn].get(name, 0.0))
        Pfix[r.bus] -= p
        gens_rows.append([bi[r.bus], p, 0.0, 0.0, 0.0, float(n.buses.v_mag_pu_set[r.bus]), base_mva, 1, max(r.p_nom, p), 0.0] + [0] * 11)
        gtype.append("HY")
    # shunts : fixes (équivalents de Kron) + commutables dans leur état de l'heure, en MVAr à 1 pu
    Bs = pd.Series(0.0, index=buses)
    sh = n.shunt_impedances
    for s, r in sh.iterrows():
        if r.bus in bi and s not in eng.sw:
            Bs[r.bus] += r.b * vn[r.bus] ** 2
    kk = K.loc[sn].reindex(buses).fillna(0.0).values; cc = C.loc[sn].reindex(buses).fillna(0.0).values
    Bs += kk * eng.b_react + cc * eng.b_cap
    # types de bus : 3 slack, 2 si un régulateur y est en PV dans la solution finale, 1 sinon
    st = S.loc[sn]
    pv_buses = {n.generators.at[gname, "bus"] for gname in eng.regs if int(st[gname]) == 0}
    Vsol = W["V"].loc[sn, buses]
    bus = []
    for b in buses:
        t = 3 if b == eng.slack_bus else (2 if b in pv_buses else 1)
        bus.append([bi[b], t, -Pfix[b], -Qfix[b], 0.0, Bs[b], 1, abs(Vsol[b]), np.degrees(np.angle(Vsol[b])), vn[b], 1, 1.1, 0.9])
    # générateurs des régulateurs basculés en PQ : Qg fixé à la limite, Qmax = Qmin = limite
    for j, gname in enumerate(eng.regs):
        s = int(st[gname])
        if s != 0:
            lim = eng.qmax[j] if s > 0 else eng.qmin[j]
            # retrouver la ligne du générateur
            busi = bi[n.generators.at[gname, "bus"]]
            for row, gt in zip(gens_rows, gtype):
                if row[0] == busi and gt == "SC" and row[1] == 0.0:
                    row[2] = lim; break
    br = []
    for l, r in n.lines.iterrows():
        if r.bus0 in bi and r.bus1 in bi:
            zb = vn[r.bus0] ** 2 / base_mva
            br.append([bi[r.bus0], bi[r.bus1], r.r / zb, r.x / zb, r.b * zb, r.s_nom, r.s_nom, r.s_nom, 0, 0, 1, -360, 360])
    for t, r in n.transformers.iterrows():
        if r.bus0 in bi and r.bus1 in bi:
            sc = base_mva / r.s_nom
            tap = float(r.get("tap_ratio", 1.0) or 1.0)
            br.append([bi[r.bus0], bi[r.bus1], r.r * sc, r.x * sc, 0.0, r.s_nom, r.s_nom, r.s_nom, tap, float(r.get("phase_shift", 0.0) or 0.0), 1, -360, 360])
    # régulateurs et slack d'abord : pandapower (from_ppc) ne garde que le 1er générateur d'un bus comme
    # générateur régulant, les suivants deviennent des injections fixes
    order = sorted(range(len(gens_rows)), key=lambda i: (0 if gtype[i] == "SC" or gens_rows[i][0] == bi[eng.slack_bus] and gens_rows[i][3] != 0 else 1, i))
    gens_rows = [gens_rows[i] for i in order]; gtype = [gtype[i] for i in order]
    return dict(version="2", baseMVA=base_mva, bus=np.array(bus, float), gen=np.array(gens_rows, float),
                branch=np.array(br, float), gentype=gtype, bus_names=list(buses))


def write_matpower(ppc, path, name="quebec_315kv_pointe", header=""):
    def mat(a):
        return "\n".join("\t" + "\t".join(f"{v:.8g}" for v in row) + ";" for row in a)
    with open(path, "w", encoding="utf-8") as f:
        f.write(f"function mpc = {name}\n%% {header}\nmpc.version = '2';\nmpc.baseMVA = {ppc['baseMVA']};\n\n")
        f.write("%% bus_i type Pd Qd Gs Bs area Vm Va baseKV zone Vmax Vmin\nmpc.bus = [\n" + mat(ppc["bus"]) + "\n];\n\n")
        f.write("%% bus Pg Qg Qmax Qmin Vg mBase status Pmax Pmin (+11 colonnes à 0)\nmpc.gen = [\n" + mat(ppc["gen"]) + "\n];\n\n")
        f.write("%% fbus tbus r x b rateA rateB rateC ratio angle status angmin angmax\nmpc.branch = [\n" + mat(ppc["branch"]) + "\n];\n\n")
        f.write("%% types de générateurs pour AC-CFM (SC = régulateur : centrales + SVC/CS du poste agrégés)\nmpc.gentype = {\n"
                + "\n".join(f"\t'{t}';" for t in ppc["gentype"]) + "\n};\n\n")
        f.write("%% index MATPOWER -> bus PyPSA\n%% " + " ".join(f"{i + 1}:{b}" for i, b in enumerate(ppc["bus_names"])) + "\n")


def ppc_ybus(ppc):
    """Ybus (pu base baseMVA) recalculée depuis le cas MATPOWER, pour vérifier l'export."""
    nb = len(ppc["bus"]); Y = np.zeros((nb, nb), complex)
    for f, t, r, x, b, *_rest in ppc["branch"]:
        tap = _rest[5] or 1.0
        i, j = int(f) - 1, int(t) - 1
        y = 1 / complex(r, x)
        Y[i, i] += (y + 1j * b / 2) / tap ** 2; Y[j, j] += y + 1j * b / 2
        Y[i, j] -= y / tap; Y[j, i] -= y / tap
    Y += np.diag(1j * ppc["bus"][:, 5] / ppc["baseMVA"])
    return Y


def check_pandapower(ppc, warm=False):
    """Troisième solveur : pandapower. 1) état final figé (bus PV/PQ de la solution, régulateurs basculés
    en injections Q fixes), départ DC ; 2) relance à chaud avec enforce_q_lims=True : pandapower vérifie
    lui-même les limites Q. Retourne les tensions et le nombre de générateurs PV hors limites."""
    import pandapower as pp
    from pandapower.converter.pypower import from_ppc
    net = from_ppc({k: ppc[k] for k in ("version", "baseMVA", "bus", "gen", "branch")}, f_hz=60, validate_conversion=False)
    out = {}
    try:
        if warm:   # départ depuis Vm, Va du cas (= solution du moteur)
            pp.runpp(net, enforce_q_lims=False, init="auto", init_vm_pu=ppc["bus"][:, 7], init_va_degree=ppc["bus"][:, 8],
                     max_iteration=50, tolerance_mva=1e-6)
        else:
            pp.runpp(net, enforce_q_lims=False, init="dc", max_iteration=50, tolerance_mva=1e-6)
        out["converge_etat_fige"] = bool(net.converged)
    except Exception:
        out["converge_etat_fige"] = False
        return out
    G = net.gen.assign(q=net.res_gen.q_mvar.values).groupby("bus")[["q", "max_q_mvar", "min_q_mvar"]].sum()
    out["bus_PV_hors_limites"] = int(((G.q > G.max_q_mvar + 1) | (G.q < G.min_q_mvar - 1)).sum())
    v1 = net.res_bus.vm_pu.values.copy()
    try:
        pp.runpp(net, enforce_q_lims=True, init="results", max_iteration=50, tolerance_mva=1e-6)
        out["converge_limites_Q_pandapower"] = bool(net.converged)
    except Exception:
        out["converge_limites_Q_pandapower"] = False
    out["v"] = v1
    out["pertes_MW"] = float(net.res_line.pl_mw.sum() + net.res_trafo.pl_mw.sum())
    return out
