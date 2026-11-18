"""Partie B : une hypothèse changée à la fois, semaine complète + marge à la pointe.
Usage : python b_variantes.py i n   (variantes i modulo n) ; écrit b_variantes_part<i>.csv
Aucune modification de reseau315.py : les changements passent par les arguments de R.build ou par des
remplacements de fonctions en mémoire (« monkeypatch »), décrits pour chaque variante."""
import sys, time
from b_outils import *

full = R.load(os.path.join(R.CURRENT, "elec_full.nc"))
T0 = R.mercator_table(full)
REF = rd("semaine_resume.csv")
LOW = [pd.Timestamp("2022-01-03 17:00"), pd.Timestamp("2022-01-03 18:00"), pd.Timestamp("2022-01-04 10:00")]
_lopf0, _regions0 = R.lopf, dict(R.KRON_REGIONS)

def T_hybride():
    """Table où `factor` (utilisé seulement par kron_equivalent) suit la règle hybride de fix_lengths."""
    T = T0.copy(); sm = float(T.loc[T.reliable, "sinuosity"].median())
    T["factor"] = np.where(T.reliable, T.true_km, sm * T.gc_km) / T.length_stored
    return T

def lopf_bruit(seed, eps=1e-3):
    def f(n, solver="highs", log=None):
        rng = np.random.default_rng(seed)
        for c in (n.generators, n.storage_units):
            c["marginal_cost"] = c["marginal_cost"] * (1 + eps * rng.uniform(-1, 1, len(c))) + eps * rng.uniform(0, 1, len(c))
        return _lopf0(n, solver=solver, log=log)
    return f

def lopf_charges(facteur):
    """Charges × facteur horaire avant le LOPF (hypothèse : la demande HQ contient déjà les pertes)."""
    def f(n, solver="highs", log=None):
        n.loads_t.p_set = n.loads_t.p_set.mul(facteur.reindex(n.snapshots).values, axis=0)
        return _lopf0(n, solver=solver, log=log)
    return f

def participation(eng, sn):
    n = eng.n; gen = np.zeros(eng.nb)
    su = n.storage_units; sp = n.storage_units_t.p_set.loc[sn].reindex(su.index).fillna(0.0)
    for nm, p in sp.items():
        if su.at[nm, "bus"] in eng.ix and p > 0: gen[eng.ix[su.at[nm, "bus"]]] += p
    g = n.generators; gp = n.generators_t.p_set.loc[sn].reindex(g.index).fillna(0.0)
    for nm, p in gp.items():
        if g.at[nm, "carrier"] in ("hydro", "ror", "OCGT") and g.at[nm, "bus"] in eng.ix and p > 0: gen[eng.ix[g.at[nm, "bus"]]] += p
    return gen / gen.sum()

def repartir_pertes(eng, passes=2):
    """Pertes AC réparties sur les centrales hydrauliques en marche (au prorata de leur production) au lieu
    d'être toutes produites par le slack. On itère : semaine -> p_slack -> répartition."""
    for it in range(passes):
        Tt = R.solve_week_reserve(eng)[0]
        for sn in Tt.index[Tt.newton_ok.astype(bool)]:
            part = participation(eng, sn); extra = Tt.p_slack[sn] * (1 - 0)   # tout ce que le slack ajoute
            add = extra * part; add[eng.sl] = 0.0     # la part du bus slack reste au slack
            eng.P.loc[sn] = eng.P.loc[sn].values + add
    return eng

def q_engage(eng, sn):
    """Limites Q de l'heure si seules les machines en marche (P > 1 MW) fournissent du réactif."""
    n = eng.n; lo = {}; hi = {}
    su = n.storage_units; sp = n.storage_units_t.p_set.loc[sn].reindex(su.index).fillna(0.0)
    g = n.generators; gp = n.generators_t.p_set.loc[sn].reindex(g.index).fillna(0.0)
    def add(b, a, c): lo[b] = lo.get(b, 0.0) + a; hi[b] = hi.get(b, 0.0) + c
    tot_inst = tot_on = 0.0
    for nm, r in su.iterrows():
        tot_inst += r.p_nom
        if sp[nm] > 1: add(r.bus, R.Q_SYNC[0] * r.p_nom, R.Q_SYNC[1] * r.p_nom); tot_on += r.p_nom
    for nm, r in g.iterrows():
        if r.carrier in R.SYNC:
            tot_inst += r.p_nom
            if gp[nm] > 1: add(r.bus, R.Q_SYNC[0] * r.p_nom, R.Q_SYNC[1] * r.p_nom); tot_on += r.p_nom
        elif r.carrier in R.RENEW and gp[nm] > 1: add(r.bus, R.Q_RENEW[0] * r.p_nom, R.Q_RENEW[1] * r.p_nom)
    for _, d in R.reactive_devices_table().iterrows(): add(d.bus, d.nombre * d.qmin_unit, d.nombre * d.qmax_unit)
    qmin = np.array([lo.get(eng.buses[b], 0.0) for b in eng.reg_bus]); qmax = np.array([hi.get(eng.buses[b], 0.0) for b in eng.reg_bus])
    return qmin, qmax, tot_inst, tot_on

VARIANTES = [
    ("Cas final (référence)", {}, None),
    ("B2 Kron : longueurs des lignes éliminées par la règle hybride", dict(T="hyb"), None),
    ("B2 Kron étendu : Saguenay + tout le sud (lat 44,9-47,6)", dict(regions={"Saguenay": R.SAGUENAY_BOX, "Sud": dict(lat=(44.9, 47.6), lon=(-77.5, -69.0))}), None),
    ("B2 Kron étendu : tout le Québec, une seule région", dict(regions={"Québec": dict(lat=(44.0, 63.0), lon=(-80.0, -56.0))}), None),
    ("B2 Kron étendu : tout le Québec + règle hybride", dict(regions={"Québec": dict(lat=(44.0, 63.0), lon=(-80.0, -56.0))}, T="hyb"), None),
    ("B4 LOPF : coûts perturbés de 0,1 % (graine 1)", dict(lopf=lopf_bruit(1)), None),
    ("B4 LOPF : coûts perturbés de 0,1 % (graine 2)", dict(lopf=lopf_bruit(2)), None),
    ("B4 LOPF : coûts perturbés de 0,1 % (graine 3)", dict(lopf=lopf_bruit(3)), None),
    ("B4 LOPF : coûts perturbés de 0,1 % (graine 4)", dict(lopf=lopf_bruit(4)), None),
    ("B5 pertes réparties sur les centrales en marche (au lieu de LG-2 seul)", {}, "pertes"),
    ("B5 charges réduites des pertes (demande HQ = charges + pertes), pertes au slack", dict(lopf="charges"), None),
    ("B5 charges réduites des pertes ET pertes réparties", dict(lopf="charges"), "pertes"),
    ("A1b sans la correction des circuits E3, pertes au slack", dict(build=dict(circuits=False)), None),
    ("B5 + A1b : charges réduites, pertes réparties, sans E3", dict(lopf="charges", build=dict(circuits=False)), "pertes"),
    ("B5 + B2 : charges réduites, pertes réparties, Kron tout le Québec + règle hybride", dict(lopf="charges", regions={"Québec": dict(lat=(44.0, 63.0), lon=(-80.0, -56.0))}, T="hyb"), "pertes"),
]

# ordre de calcul : les plus importantes d'abord (le temps est compté)
ORDRE = [0, 9, 10, 11, 1, 2, 12, 4, 5, 13, 6, 14, 7, 3, 8]
VARIANTES = [VARIANTES[i] for i in ORDRE]

def run(label, opt, post):
    t = time.time()
    R.lopf = _lopf0; R.KRON_REGIONS = dict(_regions0)
    if opt.get("regions"): R.KRON_REGIONS = opt["regions"]
    if opt.get("lopf") == "charges":
        fac = 1 - (REF.p_slack / eng_ref_load).clip(0, 0.1); R.lopf = lopf_charges(fac)
    elif opt.get("lopf"): R.lopf = opt["lopf"]
    T = T_hybride() if opt.get("T") == "hyb" else T0
    with R.quiet():
        n, info = R.build(T=T, full=None, **opt.get("build", {}))
    R.lopf = _lopf0; R.KRON_REGIONS = dict(_regions0)
    eng = R.ACEngine(n, info["ac"]["regulators"])
    if post == "pertes": repartir_pertes(eng)
    Tt, V, Q, K, C, S = R.solve_week_reserve(eng)
    nw = Tt.newton_ok.astype(bool); ok = Tt.converged.astype(bool)
    hard = Tt.in_hard.fillna(False).astype(bool) & ok; band = Tt.in_band.fillna(False).astype(bool) & ok
    load = n.loads_t.p_set.sum(axis=1)
    row = dict(variante=label, lignes=len(n.lines), transfos=len(n.transformers), delestage_pct=round(info["lopf"]["shed_pct"], 4),
               newton=int(nw.sum()), conv_limQ=int(ok.sum()), dans_090_110=int(hard.sum()), dans_095_105=int(band.sum()),
               v_min=round(float(Tt.v_min[ok].min()), 3) if ok.any() else np.nan, v_max=round(float(Tt.v_max[ok].max()), 3) if ok.any() else np.nan,
               charge_pointe_GW=round(load[PEAK] / 1e3, 2), pslack_pointe=round(float(Tt.p_slack.get(PEAK, np.nan))), pslack_max=round(float(Tt.p_slack.max())),
               react_max_GVAr=round(float(Tt.q_react_MVAr.max()) / 1e3, 1), react_pointe_GVAr=round(float(Tt.q_react_MVAr.get(PEAK, np.nan)) / 1e3, 1))
    for h in LOW: row[f"vmin_{h:%d_%Hh}"] = round(float(Tt.v_min.get(h, np.nan)), 3)
    row["n_pq_pointe"] = int(Tt.n_pq.get(PEAK, -1)) if ok.get(PEAK, False) else -1
    if nw.get(PEAK, False):
        d1 = R.load_margin(eng, PEAK, mode="remanoeuvre"); d2 = R.load_margin(eng, PEAK, mode="remanoeuvre", compl_tol=9.0)
        row["marge_stricte"] = float(d1.lam.max()) if len(d1) else np.nan
        row["marge_limQ"] = float(d2.lam.max()) if len(d2) else np.nan
        row["marge_limQ_V090"] = float(d2.lam[d2.v_min >= 0.90].max()) if len(d2) and (d2.v_min >= 0.90).any() else np.nan
    row["s"] = round(time.time() - t)
    return row

if __name__ == "__main__":
    i, k = int(sys.argv[1]), int(sys.argv[2])
    src = R.load(os.path.join(R.CURRENT, "elec_solved.nc")); eng_ref_load = src.loads_t.p_set.sum(axis=1).reindex(REF.index)
    out = f"b_variantes_part{i}.csv"; done = pd.read_csv(out).variante.tolist() if os.path.exists(out) else []
    only = sys.argv[3:] 
    for j, (label, opt, post) in enumerate(VARIANTES):
        if j % k != i or label in done: continue
        try: row = run(label, opt, post)
        except Exception as e:
            import traceback; row = dict(variante=label, erreur=repr(e)[:300]); traceback.print_exc()
        print(row, flush=True)
        pd.DataFrame([row]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
