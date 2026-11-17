"""
inspect_315kv.py

Inspection, de fond en comble, du réseau 315 kV « convergé » d'Etienne, reproduit par
reproduce_etienne.py (results/elec_315kv_shunt_pf.nc, 162/168 convergés en PyPSA 1.3).

Aucune nouvelle hypothèse n'est ajoutée au modèle : on lit le réseau final et l'AC PF
déjà calculé. La seule simulation est un DC power flow (n.lpf) sur une copie du réseau
final. Les sorties AC des heures divergées ne veulent rien dire, et ce DC PF sert à
estimer, pour ces heures-là, l'angle de transport et le rapport P / Pmax de chaque ligne,
avec Pmax = V^2 / x, la limite de transfert d'une ligne sans perte à 1 pu aux deux bouts.

Sorties (results/inspection/) :
    etat_des_lieux_315kv.md   rapport lisible, à donner au projet Claude
    *.csv                     tables complètes derrière chaque section

Usage (racine du dépôt) :
    pypsa-env/Scripts/python.exe Jean_tests/convergence_315kv/inspect_315kv.py
"""
import logging
import os
import re
import warnings

import numpy as np
import pandas as pd
import pypsa

warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
RES = os.path.join(HERE, "results")
OUT = os.path.join(RES, "inspection")
F_FINAL = os.path.join(RES, "elec_315kv_shunt_pf.nc")
F_BASE = os.path.join(REPO, "network", "networks_current", "elec_solved.nc")
F_DIVERGED = os.path.join(RES, "reproduction_etienne_heures_divergees.txt")

BAND = (0.90, 1.10)     # critère de Jean
BAND_OP = (0.95, 1.05)  # bande d'exploitation usuelle
Q_CAP = np.tan(np.arccos(0.9))  # |Q| <= 0,484 p_nom : capacité générique, cohérente avec run_pf.py (fp 0,9)
DELTA_WARN = 30.0       # degrés : angle de transport au-delà duquel une ligne est considérée « tendue »


# --------------------------------------------------------------------------- utilitaires
def md(df: pd.DataFrame, floatfmt: str = "{:.3g}", index: bool = True, max_rows: int = 25) -> str:
    """Tableau markdown sans dépendance (tabulate n'est pas installé)."""
    d = df.head(max_rows).copy()
    if index:
        d = d.reset_index()
    fmt = lambda x: floatfmt.format(x) if isinstance(x, (float, np.floating)) and np.isfinite(x) else str(x)  # noqa: E731
    lines = ["| " + " | ".join(map(str, d.columns)) + " |", "|" + "---|" * len(d.columns)]
    lines += ["| " + " | ".join(fmt(v) for v in row) + " |" for row in d.itertuples(index=False)]
    if len(df) > max_rows:
        lines.append(f"\n*({len(df) - max_rows} lignes de plus dans le CSV)*")
    return "\n".join(lines)


def save(df: pd.DataFrame, name: str) -> str:
    df.to_csv(os.path.join(OUT, name))
    return name


def diverged_hours(n: pypsa.Network) -> pd.DatetimeIndex:
    with open(F_DIVERGED) as fh:
        line = [l for l in fh if l.startswith("3_shunts")][0]
    stamps = re.findall(r"\d{4}-\d\d-\d\d \d\d:\d\d:\d\d", line)
    return pd.DatetimeIndex(stamps).intersection(n.snapshots)


def bus_degree(n: pypsa.Network) -> pd.Series:
    ends = pd.concat([n.lines.bus0, n.lines.bus1, n.transformers.bus0, n.transformers.bus1])
    return ends.value_counts().reindex(n.buses.index, fill_value=0)


# --------------------------------------------------------------------------- main
def main():
    os.makedirs(OUT, exist_ok=True)
    n = pypsa.Network(F_FINAL)
    base = pypsa.Network(F_BASE)
    S = n.snapshots
    bad = diverged_hours(n)
    ok = S.difference(bad)
    ac = n.buses.index[n.buses.carrier == "AC"]
    rep = []  # morceaux du rapport
    add = rep.append

    load_p = n.loads_t.p_set.T.groupby(n.loads.bus).sum().T  # MW par bus
    load_q = n.loads_t.q_set.T.groupby(n.loads.bus).sum().T
    tot_load = load_p.sum(axis=1)
    real = n.generators.carrier != "load_shedding"
    deg = bus_degree(n)

    add("# État des lieux — réseau 315 kV « convergé » d'Etienne\n")
    add(f"Fichier inspecté : `{os.path.relpath(F_FINAL, REPO)}` (reproduction de la séquence d'Etienne, "
        f"commit 2ef0c21d, PyPSA {pypsa.__version__}). Référence avant corrections : "
        f"`network/networks_current/elec_solved.nc`. Généré par `inspect_315kv.py`.\n")
    add("Chaîne de traitement : `elec_solved.nc` → `correct_sending_end_circuits.py` → "
        "`apply_series_compensation.py` (70 %, 15 lignes de Nathan + 4 ajoutées par Etienne) → "
        "`add_shunt_impedance.py` (réacteurs ×1,0, condensateurs ×0,7) → AC PF (`run_pf.py`, "
        "`--pv-min-capacity 0 --pv-max-load-ratio 1`).\n")

    # ------------------------------------------------------------------ 0. synthèse chiffrée
    v = n.buses_t.v_mag_pu.loc[ok, ac]
    in_band = ((v >= BAND[0]) & (v <= BAND[1])).all(axis=1)
    in_op = ((v >= BAND_OP[0]) & (v <= BAND_OP[1])).all(axis=1)
    add("## 0. Synthèse\n")
    add(f"- Snapshots : {len(S)} ; Newton converge sur **{len(ok)}** ; divergés : {len(bad)} "
        f"({', '.join(t.strftime('%m-%d %Hh') for t in bad)}).")
    add(f"- Heures convergées avec tous les bus AC dans {BAND[0]}-{BAND[1]} pu : **{int(in_band.sum())}** ; "
        f"dans {BAND_OP[0]}-{BAND_OP[1]} pu : **{int(in_op.sum())}**.")
    add(f"- Tension sur les heures convergées : min {v.min().min():.3f} pu, max {v.max().max():.3f} pu ; "
        f"bus AC qui sortent au moins une fois de {BAND[0]}-{BAND[1]} : "
        f"{int(((v < BAND[0]) | (v > BAND[1])).any().sum())} / {len(ac)}.")
    add(f"- Demande : moyenne {tot_load.mean():.0f} MW, pointe {tot_load.max():.0f} MW "
        f"({tot_load.idxmax():%m-%d %Hh}) ; moyenne des heures divergées {tot_load[bad].mean():.0f} MW, "
        f"des heures convergées {tot_load[ok].mean():.0f} MW.\n")

    # ------------------------------------------------------------------ 1. inventaire
    add("## 1. Inventaire du réseau\n")
    buses = n.buses.groupby(["carrier", "v_nom"]).size().rename("bus").to_frame()
    add("**Bus**\n\n" + md(buses) + "\n")
    lines = n.lines.groupby("v_nom").agg(lignes=("x", "size"), km=("length", "sum"),
                                         km_max=("length", "max"), circuits=("num_parallel", "sum"))
    add("**Lignes AC**\n\n" + md(lines, "{:.0f}") + "\n")
    add(f"**Transformateurs** : {len(n.transformers)}, s_nom = {n.transformers.s_nom.unique().tolist()} MVA "
        f"(valeur plate de `run_lopf_main_island.py`), x = {n.transformers.x.unique().round(3).tolist()} "
        f"(pu sur s_nom), r = x/30 (générique, run_pf).\n")
    add("**Liens (HVDC / B2B)**\n\n" + md(n.links[["bus0", "bus1", "carrier", "p_nom"]], "{:.0f}") + "\n")
    gen = n.generators.assign(real=real)
    p_mean = n.generators_t.p.loc[ok].mean()
    gtab = gen.groupby(["carrier", "control"]).agg(nb=("p_nom", "size"), p_nom_MW=("p_nom", "sum"))
    gtab["p_moy_conv_MW"] = p_mean.groupby([gen.carrier, gen.control]).sum()
    add("**Générateurs** (contrôle après `prepare_for_ac_pf`)\n\n" + md(gtab, "{:.0f}") + "\n")
    su = n.storage_units.groupby("carrier").agg(nb=("p_nom", "size"), p_nom_MW=("p_nom", "sum"))
    su["p_moy_conv_MW"] = n.storage_units_t.p.loc[ok].mean().groupby(n.storage_units.carrier).sum()
    add("**Stockage (réservoirs hydro)**\n\n" + md(su, "{:.0f}") + "\n")
    add(f"**Charges** : {len(n.loads)} sur {n.loads.bus.nunique()} bus ; P moyen {tot_load.mean():.0f} MW, "
        f"Q = P·tan(acos 0,95) (générique, run_pf).\n")
    sl = n.generators[n.generators.control == "Slack"][["bus", "carrier", "p_nom"]]
    add("**Slack par sous-réseau**\n\n" + md(sl, "{:.0f}") + "\n")

    # ------------------------------------------------------------------ 2. ce qui a été modifié
    add("## 2. Ce que la chaîne a modifié par rapport à `elec_solved.nc`\n")
    common = n.lines.index.intersection(base.lines.index)
    cmp = pd.DataFrame({
        "bus0": n.lines.loc[common, "bus0"], "bus1": n.lines.loc[common, "bus1"],
        "v_nom": n.lines.loc[common, "v_nom"], "km": n.lines.loc[common, "length"],
        "circ_avant": base.lines.loc[common, "num_parallel"], "circ_apres": n.lines.loc[common, "num_parallel"],
        "x_avant": base.lines.loc[common, "x"], "x_apres": n.lines.loc[common, "x"],
        "b_avant": base.lines.loc[common, "b"], "b_apres": n.lines.loc[common, "b"],
        "s_nom_avant": base.lines.loc[common, "s_nom"], "s_nom_apres": n.lines.loc[common, "s_nom"],
    })
    cmp["x_ratio"] = cmp.x_apres / cmp.x_avant
    changed = cmp[(cmp.x_ratio.sub(1).abs() > 1e-6) | (cmp.circ_avant != cmp.circ_apres)
                  | (cmp.b_avant / cmp.b_apres).sub(1).abs().gt(1e-6)].sort_values("x_ratio")
    save(changed, "lignes_modifiees.csv")
    add(f"{len(changed)} lignes modifiées (circuits, x ou b). Rapport x_après / x_avant : "
        f"0,3 = compensation série seule (70 %) ; 1/3 = passage de 1 à 3 circuits ; ≈ 0,1 = les deux.\n")
    add(md(changed[["bus0", "bus1", "v_nom", "km", "circ_avant", "circ_apres", "x_ratio",
                    "s_nom_avant", "s_nom_apres"]], "{:.3g}", max_rows=40) + "\n")
    sh = n.shunt_impedances.copy()
    sh["Q_MVAr_1pu"] = sh.b * sh.v_nom ** 2  # + = condensateur (fournit), − = réacteur (absorbe)
    vsh = n.buses_t.v_mag_pu.loc[ok, sh.bus]
    sh["Q_moy_MVAr"] = (vsh.values ** 2 * sh.Q_MVAr_1pu.values).mean(axis=0)
    sh["V_min"] = vsh.min().values
    sh["V_max"] = vsh.max().values
    save(sh, "shunts.csv")
    add(f"**Shunts fixes ajoutés** : {int((sh.b > 0).sum())} condensateurs "
        f"(+{sh.Q_MVAr_1pu[sh.b > 0].sum():.0f} MVAr à 1 pu), {int((sh.b < 0).sum())} réacteurs "
        f"({sh.Q_MVAr_1pu[sh.b < 0].sum():.0f} MVAr). Le dimensionnement a ignoré 120 bus qui sortent de "
        f"la bande dans les deux sens au cours de la semaine (un shunt fixe ne peut corriger qu'un sens).\n")
    add(md(sh[["bus", "v_nom", "Q_MVAr_1pu", "Q_moy_MVAr", "V_min", "V_max"]], "{:.3g}") + "\n")

    # ------------------------------------------------------------------ 3. tensions
    add("## 3. Tensions sur les heures convergées\n")
    vall = n.buses_t.v_mag_pu.loc[ok, ac]
    bt = pd.DataFrame({
        "v_nom": n.buses.loc[ac, "v_nom"], "controle": n.buses.loc[ac, "control"],
        "V_min": vall.min(), "V_moy": vall.mean(), "V_max": vall.max(),
        "h_sous_090": (vall < BAND[0]).sum(), "h_sur_110": (vall > BAND[1]).sum(),
        "degre": deg.reindex(ac), "charge_moy_MW": load_p.mean().reindex(ac).fillna(0),
        "gen_p_nom_MW": n.generators[real].groupby("bus").p_nom.sum().reindex(ac).fillna(0),
        "Q_fixe_genPQ_MVAr": n.generators[real & (n.generators.control == "PQ")]
                              .groupby("bus").q_set.sum().reindex(ac).fillna(0),
        "shunt_MVAr_1pu": sh.groupby("bus").Q_MVAr_1pu.sum().reindex(ac).fillna(0),
    })
    save(bt, "tensions_par_bus.csv")
    by_level = bt.groupby("v_nom").agg(bus=("V_min", "size"), V_min=("V_min", "min"), V_max=("V_max", "max"),
                                       bus_hors_bande=("h_sur_110", lambda s: int(((s > 0) | (bt.loc[s.index, "h_sous_090"] > 0)).sum())))
    add("**Par niveau de tension**\n\n" + md(by_level, "{:.3f}") + "\n")
    cols = ["v_nom", "controle", "V_min", "V_max", "h_sur_110", "h_sous_090", "degre",
            "charge_moy_MW", "gen_p_nom_MW", "Q_fixe_genPQ_MVAr", "shunt_MVAr_1pu"]
    add("**Bus les plus hauts** (V max)\n\n" + md(bt.sort_values("V_max", ascending=False)[cols], max_rows=20) + "\n")
    add("**Bus les plus bas** (V min)\n\n" + md(bt.sort_values("V_min")[cols], max_rows=15) + "\n")
    per_hour = pd.DataFrame({"charge_MW": tot_load[ok], "V_min": vall.min(axis=1), "V_max": vall.max(axis=1),
                             "bus_sur_110": (vall > BAND[1]).sum(axis=1), "bus_sous_090": (vall < BAND[0]).sum(axis=1)})
    save(per_hour, "tensions_par_heure.csv")
    c = per_hour[["charge_MW", "V_max", "V_min", "bus_sur_110"]].corr().loc["charge_MW"]
    add(f"Par heure : le nombre de bus > {BAND[1]} pu varie de {per_hour.bus_sur_110.min()} à "
        f"{per_hour.bus_sur_110.max()} ; corrélation avec la demande : V_max {c.V_max:+.2f}, "
        f"V_min {c.V_min:+.2f}, nb de bus > 1,10 {c.bus_sur_110:+.2f} (CSV `tensions_par_heure.csv`).\n")

    # Zoom sur la chaîne qui monte le plus haut : voisinage des 5 bus les plus hauts
    top = bt.sort_values("V_max", ascending=False).index[:6]
    L = n.lines[n.lines.bus0.isin(top) | n.lines.bus1.isin(top)].copy()
    L["Q_charge_MVAr"] = L.b * L.v_nom ** 2  # réactif produit par la capacité de la ligne à 1 pu
    L["P0_moy_MW"] = n.lines_t.p0.loc[ok, L.index].mean()
    L["Q0_moy_MVAr"] = n.lines_t.q0.loc[ok, L.index].mean()
    L["SIL_MW"] = L.v_nom ** 2 / np.sqrt(L.x / L.b)  # puissance naturelle (Z_c = sqrt(x/b), ligne sans perte)
    save(L, "zoom_bus_les_plus_hauts.csv")
    add("**Zoom : lignes autour des 6 bus les plus hauts**. Une ligne qui transporte moins que sa "
        "puissance naturelle (SIL) produit plus de réactif qu'elle n'en consomme, et la tension monte "
        "en bout de ligne (effet Ferranti).\n")
    add(md(L[["bus0", "bus1", "v_nom", "length", "num_parallel", "Q_charge_MVAr", "SIL_MW", "P0_moy_MW",
              "Q0_moy_MVAr"]], "{:.3g}") + "\n")
    gz = n.generators[n.generators.bus.isin(top) & real][["bus", "carrier", "control", "p_nom", "q_set"]]
    add("Générateurs sur ces bus (q_set = Q imposé si PQ) :\n\n" + md(gz, "{:.3g}") + "\n")

    # ------------------------------------------------------------------ 4. réactif
    add("## 4. Puissance réactive\n")
    q = n.generators_t.q.loc[ok]
    g = n.generators.loc[q.columns]
    rg = g.carrier != "load_shedding"
    qt = pd.DataFrame({"bus": g.bus, "carrier": g.carrier, "controle": g.control, "p_nom": g.p_nom,
                       "P_moy": n.generators_t.p.loc[ok, q.columns].mean(),
                       "Q_min": q.min(), "Q_max": q.max()})
    qt["Q_max_sur_pnom"] = qt[["Q_min", "Q_max"]].abs().max(axis=1) / qt.p_nom.replace(0, np.nan)
    save(qt, "reactif_generateurs.csv")
    pv = qt[qt.controle.isin(["PV", "Slack"]) & rg]
    over = pv[(pv.Q_max_sur_pnom > Q_CAP) | pv.p_nom.eq(0)]
    add(f"Limite générique utilisée pour juger : |Q| ≤ {Q_CAP:.3f} × p_nom (facteur de puissance 0,9, la "
        f"même hypothèse que run_pf.py). PyPSA n'impose **aucune** limite Q en PF : les PV et le slack "
        f"fournissent ce qu'il faut.\n")
    add(f"- {len(pv)} générateurs PV/Slack ; **{len(over)}** dépassent la limite au moins une heure "
        f"(ou ont p_nom = 0).")
    pq = qt[(qt.controle == "PQ") & rg]
    add(f"- {len(pq)} générateurs PQ injectent un Q **fixe** de 0,484 × p_nom, même à P = 0 : "
        f"{pq.Q_max.sum():.0f} MVAr au total, toute la semaine.")
    sl_q = q[sl.index] if len(sl) else pd.DataFrame()
    for s in sl.index:
        add(f"- Slack `{s}` (bus {n.generators.at[s, 'bus']}, p_nom {n.generators.at[s, 'p_nom']:.0f} MW) : "
            f"P {n.generators_t.p.loc[ok, s].min():.0f} à {n.generators_t.p.loc[ok, s].max():.0f} MW, "
            f"Q {sl_q[s].min():.0f} à {sl_q[s].max():.0f} MVAr.")
    add("\n**PV / Slack les plus sollicités en réactif**\n\n"
        + md(over.sort_values("Q_max_sur_pnom", ascending=False), "{:.3g}", max_rows=20) + "\n")
    # bilan réactif du système (moyenne des heures convergées)
    q_gen = q.loc[:, rg].sum(axis=1)
    q_load = load_q.loc[ok].sum(axis=1)
    q_sh = (n.buses_t.v_mag_pu.loc[ok, sh.bus].values ** 2 * sh.Q_MVAr_1pu.values).sum(axis=1)
    q_lines = -(n.lines_t.q0.loc[ok] + n.lines_t.q1.loc[ok]).sum(axis=1)  # + = lignes nettement productrices
    q_tr = -(n.transformers_t.q0.loc[ok] + n.transformers_t.q1.loc[ok]).sum(axis=1) if not n.transformers_t.q0.empty else 0
    bil = pd.DataFrame({"Q_generateurs": q_gen, "Q_charges": -q_load, "Q_shunts": q_sh,
                        "Q_net_lignes": q_lines, "Q_net_transfos": q_tr}).mean().to_frame("MVAr_moy")
    add("**Bilan réactif moyen (heures convergées)**, + = fourni au réseau. Les lignes produisent "
        "(capacité, b) et consomment (I²x) ; leur solde net est la somme q0 + q1, changée de signe.\n\n"
        + md(bil, "{:.0f}") + "\n")

    # Production réactive des lignes par niveau : capacité b·V² comparée au transit (P vs SIL)
    rows = []
    for kv in sorted(n.lines.v_nom.unique()):
        l = n.lines[n.lines.v_nom == kv]
        v0v1 = (n.buses_t.v_mag_pu.loc[ok, l.bus0].values * n.buses_t.v_mag_pu.loc[ok, l.bus1].values)
        sil = l.v_nom ** 2 / np.sqrt(l.x / l.b)  # x compensé inclus pour les lignes compensées
        p = n.lines_t.p0.loc[ok, l.index].abs().mean()
        rows.append({"v_nom": kv, "lignes": len(l), "Q_capacite_1pu_MVAr": (l.b * l.v_nom ** 2).sum(),
                     "Q_capacite_V_reel_MVAr": ((l.b * l.v_nom ** 2).values * v0v1).sum(axis=1).mean(),
                     "P_sur_SIL_median": (p / sil).median(), "part_lignes_sous_SIL": (p < sil).mean()})
    qc = pd.DataFrame(rows).set_index("v_nom")
    save(qc, "reactif_capacite_lignes.csv")
    add("**Réactif produit par la capacité des lignes** (b·V²). Une ligne chargée sous sa puissance "
        "naturelle (SIL) produit plus de réactif qu'elle n'en consomme. Le modèle n'a pas de "
        "réactances shunt de ligne sur le 735 kV ; ce surplus doit donc être absorbé par les PV et le slack.\n\n"
        + md(qc, "{:.3g}") + "\n")

    # ------------------------------------------------------------------ 5. lignes : charge et angle
    add("## 5. Lignes : charge et angle de transport (heures convergées, AC)\n")
    ang = n.buses_t.v_ang.loc[ok]
    dl = np.degrees(ang[n.lines.bus0].values - ang[n.lines.bus1].values)
    delta = pd.DataFrame(dl, index=ok, columns=n.lines.index).abs()
    s_ac = np.sqrt(n.lines_t.p0.loc[ok] ** 2 + n.lines_t.q0.loc[ok] ** 2)
    load_pct = s_ac.div(n.lines.s_nom, axis=1) * 100
    lt = pd.DataFrame({"bus0": n.lines.bus0, "bus1": n.lines.bus1, "v_nom": n.lines.v_nom,
                       "km": n.lines.length, "circ": n.lines.num_parallel,
                       "delta_max_deg": delta.max(), "charge_max_pct": load_pct.max(),
                       "P_max_MW": n.lines_t.p0.loc[ok].abs().max()})
    save(lt, "lignes_ac.csv")
    add(f"- Charge max (S/s_nom) : {load_pct.max().max():.1f} % ; s_nom = St Clair × 3, donc gonflé "
        f"(voir CLAUDE.md) : la charge thermique n'est pas le bon indicateur ici.")
    add(f"- Lignes avec un angle > {DELTA_WARN:.0f}° au moins une heure : "
        f"{int((lt.delta_max_deg > DELTA_WARN).sum())} / {len(lt)}.\n")
    add(md(lt.sort_values("delta_max_deg", ascending=False), "{:.3g}", max_rows=15) + "\n")

    # ------------------------------------------------------------------ 6. heures divergées (DC PF)
    add("## 6. Les heures divergées, vues en DC power flow\n")
    m = n.copy()
    m.lpf(S)
    pdc = m.lines_t.p0
    angd = m.buses_t.v_ang
    ddc = pd.DataFrame(np.degrees(angd[m.lines.bus0].values - angd[m.lines.bus1].values),
                       index=S, columns=m.lines.index).abs()
    pmax = m.lines.v_nom ** 2 / m.lines.x  # MW, V = 1 pu aux deux bouts, ligne sans perte
    ratio = pdc.abs().div(pmax, axis=1)
    hrs = pd.DataFrame({"charge_MW": tot_load, "rang_demande": tot_load.rank(ascending=False).astype(int),
                        "diverge": S.isin(bad), "delta_DC_max_deg": ddc.max(axis=1),
                        "ligne_delta_max": ddc.idxmax(axis=1), "P_sur_Pmax_max": ratio.max(axis=1)})
    for lid in ["3360", "3352"]:
        if lid in m.lines.index:
            hrs[f"delta_{lid}_deg"] = ddc[lid]
            hrs[f"P_{lid}_MW"] = pdc[lid]
    save(hrs, "heures_dc.csv")
    add("Pour chaque heure : demande, rang (1 = pointe), plus grand angle DC et rapport P / Pmax le plus "
        "élevé (Pmax = V²/x). En DC, sin δ ≈ δ : un angle DC de plus de ~60-70° indique qu'une solution "
        "AC n'existe sans doute pas sans soutien de tension local.\n")
    add("**Heures divergées**\n\n" + md(hrs.loc[bad], "{:.3g}") + "\n")
    add("**Comparaison convergées / divergées (moyennes)**\n\n"
        + md(hrs.groupby("diverge")[["charge_MW", "rang_demande", "delta_DC_max_deg", "P_sur_Pmax_max"]].mean(), "{:.3g}") + "\n")
    worst = pd.DataFrame({"bus0": m.lines.bus0, "bus1": m.lines.bus1, "v_nom": m.lines.v_nom, "km": m.lines.length,
                          "circ": m.lines.num_parallel, "x_ohm": m.lines.x,
                          "delta_DC_max_div": ddc.loc[bad].max(), "delta_DC_max_conv": ddc.loc[ok].max(),
                          "P_sur_Pmax_div": ratio.loc[bad].max(), "P_sur_Pmax_conv": ratio.loc[ok].max()})
    worst = worst.sort_values("delta_DC_max_div", ascending=False)
    save(worst, "lignes_dc_heures_divergees.csv")
    add("**Lignes les plus tendues aux heures divergées** (angle DC max), comparées aux heures convergées\n\n"
        + md(worst, "{:.3g}", max_rows=15) + "\n")

    # ------------------------------------------------------------------ 6b. ponts et poches
    import networkx as nx
    G = nx.MultiGraph()
    for i, r in n.lines.iterrows():
        G.add_edge(r.bus0, r.bus1, key=i)
    for i, r in n.transformers.iterrows():
        G.add_edge(r.bus0, r.bus1, key=i)
    simple = nx.Graph(G)
    bridges = [(a, b) for a, b in nx.bridges(simple) if G.number_of_edges(a, b) == 1]
    gp = n.generators_t.p
    prow = []
    for a, b in bridges:
        lid = [k for k in G[a][b]][0]
        if lid not in n.lines.index:
            continue
        H = G.copy()
        H.remove_edge(a, b, key=lid)
        pocket = min(nx.connected_components(H), key=len)
        if len(pocket) < 2:
            continue  # bus terminal isolé : moins intéressant, ignoré ici
        g = n.generators[n.generators.bus.isin(pocket) & real]
        s = n.storage_units[n.storage_units.bus.isin(pocket)]
        ld = load_p.reindex(columns=list(pocket), fill_value=0).sum(axis=1)
        prod = gp.reindex(columns=g.index).fillna(0).sum(axis=1) + \
            n.storage_units_t.p.reindex(columns=s.index).fillna(0).sum(axis=1)
        prow.append({"ligne": lid, "bus0": a, "bus1": b, "v_nom": n.lines.at[lid, "v_nom"],
                     "km": n.lines.at[lid, "length"], "bus_poche": len(pocket),
                     "charge_moy_MW": ld.mean(), "prod_installee_MW": g.p_nom.sum() + s.p_nom.sum(),
                     "prod_moy_MW": prod.mean(), "Pmax_MW": n.lines.at[lid, "v_nom"] ** 2 / n.lines.at[lid, "x"],
                     "P_DC_max_MW": pdc[lid].abs().max(), "P_sur_Pmax_max": ratio[lid].max(),
                     "s_nom_MVA": n.lines.at[lid, "s_nom"],
                     "bus": ", ".join(sorted(pocket)[:12]) + (" …" if len(pocket) > 12 else "")})
    pk = pd.DataFrame(prow).set_index("ligne").sort_values("P_sur_Pmax_max", ascending=False)
    save(pk, "ponts_et_poches.csv")
    add("## 6b. Ponts (lignes uniques dont la perte isole une poche d'au moins 2 bus)\n")
    add(f"{len(bridges)} ponts à circuit unique au total (lignes + transfos), dont {len(pk)} isolent une "
        f"poche d'au moins 2 bus. Pmax = V²/x à 1 pu ; s_nom est la limite que le LOPF a vue (St Clair × 3), "
        f"le LOPF ne voit **pas** Pmax.\n")
    add(md(pk, "{:.3g}", max_rows=20) + "\n")

    # ------------------------------------------------------------------ 7. limites connues
    add("## 7. Hypothèses génériques actives dans ce résultat (rappel)\n")
    add("- Q des charges : fp 0,95 partout. Q des générateurs PQ : 0,484 × p_nom fixe (pas lié au dispatch).\n"
        "- Aucune limite Q sur PV/Slack (PyPSA PF). Consigne de tension des PV : "
        f"{sorted(n.buses.v_mag_pu_set.round(3).unique().tolist())[:5]} pu.\n"
        "- Classement PV : `--pv-min-capacity 0 --pv-max-load-ratio 1` (réglage 735 kV, pas 315 kV).\n"
        "- Compensation série 70 % (générique) ; 4 lignes 735 kV ajoutées par Etienne (3535, 161, 1653, 2387), "
        "choisies empiriquement sur la convergence.\n"
        "- Shunts fixes dimensionnés à l'heure la pire, réacteurs ×1,0 / condensateurs ×0,7.\n"
        "- s_nom = St Clair × 3 ; transformateurs plats à 100 000 MVA ; dispatch = LOPF du réseau "
        "**non corrigé** (elec_solved), jamais ré-optimisé après les corrections.\n")

    path = os.path.join(OUT, "etat_des_lieux_315kv.md")
    with open(path, "w", encoding="utf-8") as fh:
        fh.write("\n".join(rep))
    print(f"Rapport : {os.path.relpath(path, REPO)}")
    print(f"CSV : {sorted(f for f in os.listdir(OUT) if f.endswith('.csv'))}")


if __name__ == "__main__":
    main()
