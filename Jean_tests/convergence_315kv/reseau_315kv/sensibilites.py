"""Alternatives de conception : on refait toute la semaine (LOPF + 168 h d'AC avec limites Q) en retirant
ou en changeant UNE hypothèse du cas final à la fois, et on calcule la marge de charge à la pointe.
Usage : python sensibilites.py [i n]   (partie i sur n, pour lancer n processus en parallèle)
Écrit resultats/sensibilites_part<i>.csv ; `python sensibilites.py --fusion` les fusionne dans
resultats/sensibilites.csv, dans l'ordre de VARIANTS. Reprend là où il s'est arrêté."""
import os, sys, time, glob, numpy as np, pandas as pd
import reseau315 as R

FINAL = dict()   # = valeurs par défaut de R.build() (cas final)
NATHAN_ETIENNE_19 = ["156", "162", "163", "880", "1049", "1448", "1536", "1538", "1655", "1935", "1936", "1938",
                     "2386", "2408", "2410", "3535", "161", "1653", "2387"]
VARIANTS = [
    ("Cas final", dict()),
    ("Sans correction des longueurs (E1)", dict(lengths=False)),
    ("Sans correction de Churchill Falls (E2)", dict(churchill=False)),
    ("Sans équivalents de Kron (E2b)", dict(kron=False)),
    ("Équivalent de Kron au Saguenay seulement (E2b)", dict(kron=["Saguenay"])),
    ("Sans correction des circuits de Nathan (E3)", dict(circuits=False)),
    ("Compensation série 70 % sur les 19 lignes Nathan/Etienne (E4)", dict(corridors={"Nathan+Etienne": NATHAN_ETIENNE_19}, series_fraction=0.70)),
    ("Sans compensation série (E4)", dict(series_fraction=0.0)),
    ("Sans plafond de stabilité dans le LOPF (E5)", dict(delta_max_deg=None)),
    ("Sans SVC ni compensateurs synchrones (E6)", dict(devices=False)),
    ("Sans bancs de condensateurs (E6c)", dict(alpha_banks=0.0)),
    ("Bancs de condensateurs alpha = 0,6 (E6c)", dict(alpha_banks=0.6)),
    ("Bancs de condensateurs alpha = 1,5 (E6c)", dict(alpha_banks=1.5)),
    ("Charges à cos φ 0,97 (E6)", dict(load_pf=0.97)),
    ("Consigne de tension 1,00 pu (E6)", dict(v_set=1.0)),
    ("Consigne de tension 1,04 pu (E6)", dict(v_set=1.04)),
]
PEAK = pd.Timestamp("2022-01-03 17:00")


def run(label, kw, T=None):
    t = time.time()
    args = dict(FINAL); args.update(kw)
    with R.quiet():
        n, info = R.build(T=T, **args)
    eng = R.ACEngine(n, info["ac"]["regulators"])
    Tt = R.solve_week_reserve(eng)[0]
    newton = Tt.newton_ok.astype(bool)
    ok = Tt.converged.astype(bool)                      # Newton ET limites Q respectées
    hard = Tt.get("in_hard", pd.Series(False, index=Tt.index)).fillna(False).astype(bool) & ok
    band = Tt.get("in_band", pd.Series(False, index=Tt.index)).fillna(False).astype(bool) & ok
    d = R.load_margin(eng, PEAK, mode="remanoeuvre") if ok.get(PEAK, False) else pd.DataFrame()
    return dict(variante=label,
                lopf_delestage_pct=round(info["lopf"]["shed_pct"], 3) if "lopf" in info else np.nan,
                newton=int(newton.sum()), convergees=int(ok.sum()),
                bande_090_110=int(hard.sum()), bande_095_105=int(band.sum()),
                v_min=round(float(Tt.v_min[ok].min()), 3) if ok.any() else np.nan,
                v_max=round(float(Tt.v_max[ok].max()), 3) if ok.any() else np.nan,
                inductances_pointe_GVAr=round(float(Tt.q_react_MVAr.loc[PEAK]) / 1e3, 1) if ok.get(PEAK, False) else np.nan,
                regulateurs_PQ_moy=round(float(Tt.n_pq[ok].mean()), 1) if ok.any() else np.nan,
                marge_pointe_lambda=round(float(d.lam.max()), 3) if len(d) else np.nan,
                s=round(time.time() - t))


if __name__ == "__main__":
    os.makedirs("resultats", exist_ok=True)
    if "--fusion" in sys.argv:
        parts = pd.concat([pd.read_csv(f) for f in sorted(glob.glob("resultats/sensibilites_part*.csv"))])
        order = {v: i for i, (v, _) in enumerate(VARIANTS)}
        parts = parts.drop_duplicates("variante", keep="last").sort_values("variante", key=lambda s: s.map(order))
        parts.to_csv("resultats/sensibilites.csv", index=False)
        print(parts.to_string(index=False)); sys.exit()
    i, k = (int(sys.argv[1]), int(sys.argv[2])) if len(sys.argv) > 2 else (0, 1)
    out = os.path.join("resultats", f"sensibilites_part{i}.csv")
    done = pd.read_csv(out).variante.tolist() if os.path.exists(out) else []
    T = R.mercator_table(R.load(os.path.join(R.CURRENT, "elec_full.nc")))
    for j, (label, kw) in enumerate(VARIANTS):
        if j % k != i or label in done:
            continue
        try:
            row = run(label, kw, T=T)
        except Exception as e:
            row = dict(variante=label, erreur=repr(e)[:200])
        print(row, flush=True)
        pd.DataFrame([row]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
