"""
reproduce_etienne.py

Refait, avec un script plutôt qu'un notebook, la séquence d'Etienne (commit 2ef0c21d,
etienne_tests/315_convergence.ipynb) qui mène à 161/168 snapshots convergés en AC PF
sur le réseau 315 kV :

    0. elec_solved.nc brut                                   (Etienne : 3/168)
    1. + correct_sending_end_circuits.py                     (Etienne : 10/168)
    2. + apply_series_compensation.py 70 %, 15+4 lignes      (Etienne : 157/168)
    3. + add_shunt_impedance.py (SCALE_GRID 1.0/0.7/0.4)     (Etienne : 161/168)

Les scripts de Nathan (dans leur version modifiée par Etienne) sont appelés tels quels
par subprocess ; rien n'est écrit dans network/networks_current/, tout va dans results/.

L'AC PF est refait ici, dans le même ordre que run_pf.main() (consignes LOPF -> p_set,
prepare_for_ac_pf avec --pv-min-capacity 0 --pv-max-load-ratio 1 comme Etienne, seed DC,
Newton avec use_seed=True), pour récupérer le drapeau de convergence de Newton par
snapshot. Le fichier *_pf.nc écrit est le même que celui de run_pf.py, et c'est lui que
l'étape shunts prend en entrée, comme chez Etienne.

Deux chiffres par étape (règle du projet) :
  - convergés : Newton converge sur tous les sous-réseaux du snapshot ;
  - dans la bande : convergés ET 0,90 <= V <= 1,10 pu sur tous les bus AC.
Le contrôle « réactif plausible » n'est pas encore fait ici (pas de limites Q dans PyPSA).

Usage (racine du dépôt) :
    pypsa-env/Scripts/python.exe Jean_tests/convergence_315kv/reproduce_etienne.py
"""
import os
import subprocess
import sys

import pandas as pd
import pypsa

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
NETWORK_DIR = os.path.join(REPO, "network")
RESULTS = os.path.join(HERE, "results")
SOURCE = os.path.join(NETWORK_DIR, "networks_current", "elec_solved.nc")

sys.path.insert(0, NETWORK_DIR)
import run_pf  # noqa: E402

PV_MIN_CAPACITY = 0.0   # mêmes réglages qu'Etienne
PV_MAX_LOAD_RATIO = 1.0
V_LO, V_HI = 0.90, 1.10  # bande du critère « convergé » de Jean
ETIENNE = {"0_brut": 3, "1_circuits": 10, "2_compensation": 157, "3_shunts": 161}


def run_script(name: str, *args: str) -> None:
    cmd = [sys.executable, os.path.join(NETWORK_DIR, name), *args]
    print("\n$", " ".join(os.path.relpath(c, REPO) if os.path.isabs(c) else c for c in cmd), flush=True)
    subprocess.run(cmd, check=True, cwd=REPO)


def ac_pf(path: str, out_pf: str) -> pd.Series:
    """AC PF identique à run_pf.main(); renvoie le drapeau de convergence par snapshot."""
    n = pypsa.Network(path)
    run_pf.use_lopf_dispatch_as_setpoints(n)
    n.determine_network_topology()
    run_pf.prepare_for_ac_pf(n, pv_min_capacity=PV_MIN_CAPACITY, pv_max_load_ratio=PV_MAX_LOAD_RATIO)
    n.lpf(n.snapshots)
    res = n.pf(n.snapshots, use_seed=True)
    n.export_to_netcdf(out_pf)
    conv = res.converged.all(axis=1)

    ac_buses = n.buses.index[n.buses.carrier == "AC"]
    v = n.buses_t.v_mag_pu.loc[conv, ac_buses]
    in_band = ((v >= V_LO) & (v <= V_HI)).all(axis=1)
    conv.attrs.update(
        in_band=int(in_band.sum()),
        v_min=float(v.min().min()) if len(v) else float("nan"),
        v_max=float(v.max().max()) if len(v) else float("nan"),
        worst_low=v.min().nsmallest(3).round(3).to_dict() if len(v) else {},
        worst_high=v.max().nlargest(3).round(3).to_dict() if len(v) else {},
        n_ac=len(ac_buses),
    )
    return conv


def main():
    os.makedirs(RESULTS, exist_ok=True)
    f_base = os.path.join(RESULTS, "elec_315kv_base.nc")
    f_cc = os.path.join(RESULTS, "elec_315kv.nc")
    f_shunt = os.path.join(RESULTS, "elec_315kv_shunt.nc")
    pf = lambda p: os.path.splitext(p)[0] + "_pf.nc"  # noqa: E731

    rows, diverged = [], {}

    def record(step, path, conv):
        a = conv.attrs
        rows.append({"etape": step, "converges": int(conv.sum()), "dans_bande_090_110": a["in_band"],
                     "snapshots": len(conv), "etienne": ETIENNE[step],
                     "v_min_converges": round(a["v_min"], 3), "v_max_converges": round(a["v_max"], 3),
                     "bus_v_min": a["worst_low"], "bus_v_max": a["worst_high"],
                     "fichier": os.path.relpath(path, REPO), "pypsa": pypsa.__version__})
        diverged[step] = [str(t) for t in conv.index[~conv]]
        print(f"\n>>> {step}: {int(conv.sum())}/{len(conv)} convergés, {a['in_band']} dans la bande "
              f"{V_LO}-{V_HI} pu (Etienne : {ETIENNE[step]})  V convergés {a['v_min']:.3f}-{a['v_max']:.3f}",
              flush=True)

    # 0. réseau brut (copie locale, pour ne rien écrire à côté de networks_current/)
    pypsa.Network(SOURCE).export_to_netcdf(f_base)
    record("0_brut", f_base, ac_pf(f_base, pf(f_base)))

    # 1. correction des circuits des corridors de départ
    run_script("correct_sending_end_circuits.py", "--network", SOURCE, "--output", f_cc)
    record("1_circuits", f_cc, ac_pf(f_cc, pf(f_cc)))

    # 2. compensation série 70 % (TARGET_LINES modifié par Etienne : +3535, 161, 1653, 2387)
    run_script("apply_series_compensation.py", "--network", f_cc, "--in-place")
    record("2_compensation", f_cc, ac_pf(f_cc, pf(f_cc)))

    # 3. shunts, dimensionnés à partir du *_pf.nc de l'étape 2 (comme Etienne)
    run_script("add_shunt_impedance.py", "--network", pf(f_cc), "--output", f_shunt)
    record("3_shunts", f_shunt, ac_pf(f_shunt, pf(f_shunt)))

    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(RESULTS, "reproduction_etienne.csv"), index=False)
    with open(os.path.join(RESULTS, "reproduction_etienne_heures_divergees.txt"), "w") as fh:
        for step, hours in diverged.items():
            fh.write(f"{step} ({len(hours)}): {', '.join(hours)}\n")

    print("\n=== Résumé (PyPSA %s) ===" % pypsa.__version__)
    print(table[["etape", "converges", "dans_bande_090_110", "etienne",
                 "v_min_converges", "v_max_converges"]].to_string(index=False))
    print("\nHeures divergées à la dernière étape :", ", ".join(diverged["3_shunts"]) or "aucune")


if __name__ == "__main__":
    main()
