"""
fix_line_lengths_mercator.py — longueurs des lignes recalculées (étape E1), en script autonome.

Mise à jour du 17 nov. 2026 : la version du 16 nov. au soir multipliait la longueur stockée par
(longueur géodésique / longueur Mercator) du tracé. C'est exact quand la longueur stockée est bien la
longueur Mercator du tracé (la plupart des lignes 735 kV), mais certaines lignes devenaient plus courtes
que la ligne droite entre leurs bus. La règle retenue (détails : reseau_315kv/Reseau_315kV_converge.ipynb, §3) :
  - tracé OSM fiable (extrémités à <= 5 km des bus, sinuosité tracé / vol d'oiseau entre 1,0 et 1,6) :
    longueur = longueur géodésique du tracé ;
  - sinon : longueur = vol d'oiseau × sinuosité médiane des tracés fiables (≈ 1,08) ;
  - puis length, r, x, b à l'échelle.
La logique est dans reseau_315kv/reseau315.py (fonctions mercator_table et fix_lengths) ; ce script l'appelle.

Usage (racine du dépôt) :
    python Jean_tests/convergence_315kv/fix_line_lengths_mercator.py \
        --network network/networks_current/elec_solved.nc \
        --output  Jean_tests/convergence_315kv/results/elec_solved_len.nc
"""
import argparse
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "reseau_315kv"))
import reseau315 as R  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--network", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--full", default=os.path.join(R.CURRENT, "elec_full.nc"))
    a = ap.parse_args()
    T = R.mercator_table(R.load(a.full))
    n = R.load(a.network)
    km0 = (n.lines.length * n.lines.num_parallel).groupby(n.lines.v_nom).sum()
    rep = R.fix_lengths(n, T)
    km1 = (n.lines.length * n.lines.num_parallel).groupby(n.lines.v_nom).sum()
    print(rep.statut.value_counts().to_string())
    ok = rep[rep.statut == "corrigée"]
    print(ok.groupby(["v_nom", "methode"]).facteur.agg(["size", "median"]).round(3).to_string())
    for v in km0.index:
        print(f"  {v:5.0f} kV : {km0[v]:8.0f} -> {km1[v]:8.0f} circuit-km")
    with R.quiet():
        n.export_to_netcdf(a.output)
    rep.to_csv(os.path.splitext(a.output)[0] + "_longueurs.csv")
    print("Écrit :", a.output)


if __name__ == "__main__":
    main()
