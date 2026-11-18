# -*- coding: utf-8 -*-
"""
construire.py — construit le réseau 315 kV final, résout la semaine et écrit tout dans resultats/.

    python Jean_tests/convergence_315kv/reseau_315kv/construire.py            (≈ 2-4 min)
    python Jean_tests/convergence_315kv/reseau_315kv/construire.py --valider  (+ PyPSA et pandapower, ≈ 10 min)

Fichiers écrits dans resultats/ :
  reseau_315kv_final.nc         réseau prêt pour l'AC (E1-E6), avec régulateurs et shunts commutables
  journal_etapes.csv            une ligne par étape du pipeline
  semaine_resume.csv            une ligne par heure (convergence, tensions, inductances, condensateurs...)
  tensions_pu.csv, angles_deg.csv           |V| et angle par bus et par heure (sous-réseau AC principal)
  regulateurs_Q_MVAr.csv, regulateurs_etat.csv   Q de chaque régulateur ; 0 = PV, +1/-1 = PQ en butée haute/basse
  inductances_k.csv, condensateurs_c.csv    part enclenchée par bus et par heure (0 à 1)
  longueurs_corrigees.csv, kron_branches.csv, equipements_reactifs.csv, capacites_Q_par_bus.csv
  quebec_315kv_pointe.m          cas MATPOWER de la pointe (3 janv. 2022 17 h) pour AC-CFM
  validation_pypsa.csv, validation_pandapower.csv  (avec --valider)
"""
import argparse
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import reseau315 as R  # noqa: E402

OUT = os.path.join(HERE, "resultats")
PEAK = pd.Timestamp("2022-01-03 17:00")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--valider", action="store_true", help="revérifier chaque heure avec PyPSA et pandapower")
    a = ap.parse_args()
    os.makedirs(OUT, exist_ok=True)
    t0 = time.time()
    print("Table Mercator (longueurs vraies vs stockées)...", flush=True)
    full = R.load(os.path.join(R.CURRENT, "elec_full.nc"))
    T = R.mercator_table(full)
    print("Construction E1-E6...", flush=True)
    n, info = R.build(T=T)
    for step, msg in info["log"]:
        print(f"  {step} : {msg}")
    pd.DataFrame(info["log"], columns=["etape", "resultat"]).to_csv(os.path.join(OUT, "journal_etapes.csv"), index=False)
    info["lengths"].to_csv(os.path.join(OUT, "longueurs_corrigees.csv"))
    kb = []
    for reg, (br, sh, kinfo) in info.get("kron", {}).items():
        kb.append(br.assign(region=reg, r_pu_1MVA=br.z_pu.apply(lambda z: z.real), x_pu_1MVA=br.z_pu.apply(lambda z: z.imag),
                            x_ohm_vu_du_315kV=br.z_pu.apply(lambda z: z.imag) * 315 ** 2).drop(columns="z_pu"))
    if kb:
        pd.concat(kb).to_csv(os.path.join(OUT, "kron_branches.csv"), index=False)
    R.reactive_devices_table().to_csv(os.path.join(OUT, "equipements_reactifs.csv"), index=False)
    info["ac"]["capability"].groupby(["bus", "source"])[["qmin", "qmax"]].sum().to_csv(os.path.join(OUT, "capacites_Q_par_bus.csv"))
    with R.quiet():
        n.export_to_netcdf(os.path.join(OUT, "reseau_315kv_final.nc"))
    print(f"  ({time.time() - t0:.0f} s)\nRésolution des 168 heures...", flush=True)
    eng = R.ACEngine(n, info["ac"]["regulators"])
    Tt, V, Q, K, C, S = R.solve_week_reserve(eng)
    W = dict(T=Tt, V=V, Q=Q, K=K, C=C, S=S)
    Tt.to_csv(os.path.join(OUT, "semaine_resume.csv"))
    V.abs().to_csv(os.path.join(OUT, "tensions_pu.csv"))
    pd.DataFrame(np.degrees(np.angle(V.values.astype(complex))), index=V.index, columns=V.columns).to_csv(os.path.join(OUT, "angles_deg.csv"))
    Q.to_csv(os.path.join(OUT, "regulateurs_Q_MVAr.csv")); S.to_csv(os.path.join(OUT, "regulateurs_etat.csv"))
    K.loc[:, eng.b_react != 0].to_csv(os.path.join(OUT, "inductances_k.csv"))
    C.loc[:, eng.b_cap != 0].to_csv(os.path.join(OUT, "condensateurs_c.csv"))
    ok = Tt.converged.astype(bool)
    print(f"  {int(ok.sum())}/168 convergées | 0,90-1,10 pu : {int(Tt.in_hard.fillna(False).astype(bool).sum())} | "
          f"0,95-1,05 pu : {int(Tt.in_band.fillna(False).astype(bool).sum())} | V {Tt.v_min.min():.3f}-{Tt.v_max.max():.3f} "
          f"({time.time() - t0:.0f} s)", flush=True)
    if ok.get(PEAK, False):
        ppc = R.to_ppc(n, eng, W, PEAK)
        R.write_matpower(ppc, os.path.join(OUT, "quebec_315kv_pointe.m"),
                         header="Réseau 315 kV du Québec (modèle de Nathan Guan corrigé, capstone ECSE 458), pointe du 3 janv. 2022 17 h, "
                                "état final (shunts commutés, limites Q). Généré par construire.py.")
        print("  cas MATPOWER de la pointe écrit", flush=True)
    if a.valider:
        print("Validation PyPSA (chaque heure, départ DC à froid)...", flush=True)
        D = R.validate_with_pypsa(n, eng, W, verbose=True)
        D.to_csv(os.path.join(OUT, "validation_pypsa.csv"))
        far = D.ecart_V_max > 1e-4
        print(f"  PyPSA (départ à froid) : {int(D.pypsa_converge.sum())}/{len(D)} | même solution (écart < 1e-4 pu) : {int((~far).sum())}/{len(D)}", flush=True)
        if far.any():
            print(f"  heures où PyPSA à froid trouve une autre solution : {[str(s) for s in D.index[far]]} ; au départ à chaud, "
                  f"écart max {D.loc[far, 'ecart_V_max_depart_chaud'].max():.1e} pu", flush=True)
        try:
            import pandapower  # noqa: F401
        except ImportError:
            pandapower = None
            print("pandapower non installé : validation pandapower sautée.", flush=True)
        print("Validation pandapower (chaque heure)...", flush=True)
        rows = []
        for sn in (Tt.index[ok] if pandapower is not None else []):
            ppc = R.to_ppc(n, eng, W, sn)
            r = R.check_pandapower(ppc)
            v = np.abs(V.loc[sn, eng.buses].values.astype(complex))
            row = dict(snapshot=sn, converge=r.get("converge_etat_fige"), bus_PV_hors_limites=r.get("bus_PV_hors_limites"),
                       limites_Q_pandapower_ok=r.get("converge_limites_Q_pandapower"),
                       ecart_V_max=float(np.abs(r["v"] - v).max()) if "v" in r else np.nan, pertes_MW=r.get("pertes_MW"))
            if not (row["ecart_V_max"] <= 1e-4):
                r2 = R.check_pandapower(ppc, warm=True)
                row.update(depart_chaud_converge=r2.get("converge_etat_fige"), depart_chaud_bus_PV_hors_limites=r2.get("bus_PV_hors_limites"),
                           ecart_V_max_depart_chaud=float(np.abs(r2["v"] - v).max()) if "v" in r2 else np.nan)
            rows.append(row)
        if rows:
            P = pd.DataFrame(rows).set_index("snapshot")
            P.to_csv(os.path.join(OUT, "validation_pandapower.csv"))
            print(f"  pandapower : {int(P.converge.sum())}/{len(P)} | écart |V| max {P.ecart_V_max.max():.1e} pu | "
                  f"limites Q revérifiées par pandapower : {int(P.limites_Q_pandapower_ok.sum())}/{len(P)}", flush=True)
    print(f"Terminé en {time.time() - t0:.0f} s -> {OUT}")


if __name__ == "__main__":
    main()
