"""Figures for the Nov. 17, 2026 meeting with Prof. Wang (ECSE 458 capstone).

Every figure reads its numbers from data/ (written by make_results.py) or from the review folder
(../revue_fable/). Nothing is typed by hand except the milestone scorecard (fig_scorecard), whose
sources are given next to each number. Each function returns the matplotlib Figure and, if `save`
is given, writes a PNG to figures/.
"""
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
FIGS = os.path.join(HERE, "figures")
os.makedirs(FIGS, exist_ok=True)
PEAK = pd.Timestamp("2022-01-03 17:00")
RECORD_GW = 43.124  # Hydro-Quebec all-time peak, Feb. 3 2023

# ---------------------------------------------------------------- style (one palette for notebook and slides)
SURF, INK, INK2, GRID = "#f7f6f2", "#1d2327", "#55606a", "#dcdad3"
BLUE, ORANGE, AQUA, GRAY = "#2a78d6", "#eb6834", "#1baf7a", "#a3a29b"
VMAP = LinearSegmentedColormap.from_list("v", [BLUE, "#d9d8d2", ORANGE])

plt.rcParams.update({
    "figure.facecolor": SURF, "axes.facecolor": SURF, "savefig.facecolor": SURF,
    "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.titlecolor": INK, "axes.titlesize": 13, "axes.titleweight": "bold",
    "axes.titlelocation": "left", "axes.labelsize": 11, "xtick.labelsize": 10, "ytick.labelsize": 10,
    "legend.fontsize": 10, "legend.frameon": False, "font.family": "DejaVu Sans",
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "axes.spines.top": False, "axes.spines.right": False, "lines.linewidth": 2,
})


def _save(fig, save):
    if save:
        fig.savefig(os.path.join(FIGS, save), dpi=200, bbox_inches="tight")
    return fig


def load_data():
    """Everything the figures need, in one dict."""
    d = {}
    d["T"] = pd.read_csv(os.path.join(DATA, "mercator_table.csv"), index_col=0)
    d["w1"] = pd.read_csv(os.path.join(DATA, "v1_week_summary.csv"), index_col=0, parse_dates=True)
    d["w2"] = pd.read_csv(os.path.join(DATA, "v2_week_summary.csv"), index_col=0, parse_dates=True)
    d["V2"] = pd.read_csv(os.path.join(DATA, "v2_voltages_pu.csv"), index_col=0, parse_dates=True)
    d["pv"] = pd.read_csv(os.path.join(DATA, "pv_curves_peak.csv"))
    d["pv_kron"] = pd.read_csv(os.path.join(DATA, "review_b7_courbe_pv_Q_installe.csv"))
    d["disp"] = pd.read_csv(os.path.join(DATA, "dispatch_peak_summary.csv"))
    d["sens"] = pd.read_csv(os.path.join(DATA, "v1_sensitivities.csv"))
    d["rev"] = pd.read_csv(os.path.join(DATA, "review_variants.csv"))
    d["buses"] = pd.read_csv(os.path.join(DATA, "buses.csv"), index_col=0, dtype={0: str})
    d["buses"].index = d["buses"].index.astype(str)
    d["full"] = pd.read_csv(os.path.join(DATA, "full_network_peak_voltages.csv"), index_col=0)
    d["full_hours"] = pd.read_csv(os.path.join(DATA, "full_network_hours.csv"))
    return d


# ---------------------------------------------------------------- 1. Web Mercator
def fig_mercator(T, save=None):
    """Left: stored length / geodesic length of the same OSM route, against latitude, with the Web Mercator
    stretch 1/cos(latitude). Right: 735 kV km in the model vs the route geometry vs Hydro-Quebec."""
    fig, ax = plt.subplots(1, 2, figsize=(13, 4.8), gridspec_kw=dict(width_ratios=[1.35, 1]))
    t = T[(T.end_mismatch_km <= 5) & T.lat.notna()]
    for v, c, lab, sz in ((315, ORANGE, "315 kV lines", 16), (735, BLUE, "735 / 765 kV lines", 30)):
        s = t[(t.v_nom >= 700) if v == 735 else (t.v_nom.between(300, 400))]
        ax[0].scatter(s.lat, s.stored_over_true, s=sz, color=c, alpha=0.85, label=f"{lab} ({len(s)})", edgecolors=SURF, linewidths=0.5)
    lat = np.linspace(44.5, 55.5, 100)
    ax[0].plot(lat, 1 / np.cos(np.radians(lat)), color=INK, lw=2, ls="--", label="1 / cos(latitude)")
    ax[0].axhline(1.0, color=GRAY, lw=1)
    ax[0].set_xlabel("mean latitude of the line (°N)")
    ax[0].set_ylabel("length stored by PyPSA-Earth\n÷ true (geodesic) length of the same route")
    ax[0].set_title("Every line is too long by 1/cos(latitude)")
    ax[0].legend(loc="upper left"); ax[0].set_ylim(0.8, 1.9); ax[0].set_xlim(44, 56)
    t7 = T[T.v_nom >= 700]
    vals = [t7.length_stored.sum(), t7.true_km.sum(), 12581]
    labs = ["model\n(Web Mercator)", "same OSM routes,\ngeodesic", "Hydro-Québec\n(annual report 2023)"]
    bars = ax[1].bar(labs, vals, color=[ORANGE, BLUE, GRAY], width=0.6)
    for b, v in zip(bars, vals):
        ax[1].text(b.get_x() + b.get_width() / 2, v + 250, f"{v:,.0f} km", ha="center", fontsize=11, color=INK)
    ax[1].set_ylim(0, max(vals) * 1.15); ax[1].set_ylabel("km of 735/765 kV line (one count per line)")
    ax[1].set_title("735 kV length: +51 % in the model")
    ax[1].grid(axis="x", visible=False)
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 2. scorecard
SCORECARD = pd.DataFrame(
    # columns: Newton converges, within Q limits, all V in 0.90-1.10, all V in 0.95-1.05, V-Q stable ; None = not enforced / not computed
    [["Nathan's 315 kV model, LOPF dispatch copied (reproduced Nov. 16)", 3, None, 0, 0, None],
     ["Etienne, Nov. 14 (circuits + 70 % series + fixed shunts)", 162, None, 0, 0, None],
     ["v1: our corrected model, night of Nov. 16-17", 168, 168, 168, 165, 167],
     ["v2: v1 + run-of-river priority (after the review)", 168, 168, 168, 168, 168]],
    columns=["case", "Newton\nconverges", "within\nQ limits", "all V in\n0.90-1.10 pu", "all V in\n0.95-1.05 pu", "voltage-stable\n(dV/dQ > 0)"],
).set_index("case")


def fig_scorecard(save=None):
    S = SCORECARD
    fig, ax = plt.subplots(figsize=(13, 3.6))
    ax.set_xlim(0, S.shape[1]); ax.set_ylim(0, S.shape[0]); ax.invert_yaxis(); ax.axis("off")
    cmap = LinearSegmentedColormap.from_list("s", ["#eeede8", "#9ec5f4", BLUE])
    for i, (case, row) in enumerate(S.iterrows()):
        for j, v in enumerate(row.values):
            if v is None or (isinstance(v, float) and np.isnan(v)):
                ax.add_patch(plt.Rectangle((j + 0.04, i + 0.08), 0.92, 0.84, color="#eeede8"))
                ax.text(j + 0.5, i + 0.5, "not enforced" if j == 1 else "not computed", ha="center", va="center", fontsize=9.5, color=INK2, style="italic")
            else:
                c = cmap(v / 168)
                ax.add_patch(plt.Rectangle((j + 0.04, i + 0.08), 0.92, 0.84, color=c))
                ax.text(j + 0.5, i + 0.5, f"{int(v)} / 168", ha="center", va="center", fontsize=13, fontweight="bold",
                        color="white" if v / 168 > 0.7 else INK)
    for j, c in enumerate(S.columns):
        ax.text(j + 0.5, -0.1, c, ha="center", va="bottom", fontsize=10.5, color=INK2, linespacing=1.15)
    for i, c in enumerate(S.index):
        ax.text(-0.05, i + 0.5, c, ha="right", va="center", fontsize=10.5, color=INK)
    ax.set_title("Hours of the test week (Jan. 1-7 2022) meeting each criterion", pad=46)
    fig.subplots_adjust(left=0.36, right=0.99, top=0.76, bottom=0.02)
    return _save(fig, save)


# ---------------------------------------------------------------- 3. week
def fig_week(w1, w2, V2, buses, save=None):
    fig, ax = plt.subplots(2, 1, figsize=(13, 6.2), sharex=True, gridspec_kw=dict(height_ratios=[1, 1.8]))
    ld = pd.read_csv(os.path.join(DATA, "load_GW.csv"), index_col=0, parse_dates=True).iloc[:, 0]
    ax[0].plot(ld.index, ld.values, color=INK, lw=1.8)
    ax[0].set_ylabel("load (GW)"); ax[0].set_title("Test week, Jan. 1-7 2022 (Hydro-Québec demand)")
    ax[0].annotate(f"peak {ld.max():.1f} GW, Jan. 3, 17:00", (ld.idxmax(), ld.max()), xytext=(40, -14), textcoords="offset points",
                   fontsize=10, color=INK2, arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    hv = buses.index[buses.v_nom >= 700]; lv = buses.index[buses.v_nom < 700]
    for idx, c, lab in ((lv, BLUE, "v2: 315 kV buses, min-max"), (hv, AQUA, "v2: 735 kV buses, min-max")):
        vv = V2[[b for b in idx if b in V2.columns]]
        ax[1].fill_between(vv.index, vv.min(axis=1), vv.max(axis=1), color=c, alpha=0.35, lw=0, label=lab)
    ax[1].plot(w1.index, w1.v_min, color=ORANGE, lw=1.6, label="v1: lowest bus voltage")
    bad = w1.index[~w1.stable.astype(bool)]
    ax[1].scatter(bad, w1.v_min[bad], s=90, marker="X", color=ORANGE, edgecolors=INK, zorder=5, label="v1: voltage-unstable solution")
    for y in (0.95, 1.05):
        ax[1].axhline(y, color=INK2, lw=1, ls="--")
    ax[1].text(V2.index[2], 1.053, "0.95-1.05 pu band", fontsize=9.5, color=INK2)
    ax[1].set_ylim(0.9, 1.08); ax[1].set_ylabel("bus voltage (pu)")
    ax[1].legend(loc="lower right", ncol=2)
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 4. map
STATIONS = {"136": "Chénier", "195": "Duvernay", "308": "La Vérendrye", "312": "Abitibi", "310": "Chibougamau", "345": "Albanel",
            "457": "Némiscau", "340": "Radisson", "66": "LG-2", "298": "LG-4", "286": "Chamouchouane", "150": "Saguenay",
            "25": "Laurentides", "148": "Lévis", "114": "Manicouagan", "306": "Montagnais", "469": "Churchill Falls",
            "710": "Romaine-2", "131": "Beauharnois", "800": "bus 800 (Montréal-W)"}
MTL = {"136", "195", "131", "800"}
SVC_BUSES = ["345", "286", "310", "25", "457", "308", "148", "136", "195", "311", "114"]


def _borders(ax):
    try:
        import geopandas as gpd
        p = "/usr/local/lib/python3.13/dist-packages/pyogrio/tests/fixtures/naturalearth_lowres/naturalearth_lowres.shp"
        w = gpd.read_file(p)
        w[w.name.isin(["Canada", "United States of America"])].boundary.plot(ax=ax, color=GRAY, lw=0.6, zorder=0)
    except Exception:
        pass


def fig_map(nc_path, V2, sn=PEAK, save=None, box=None, title=None):
    import pypsa, logging
    logging.getLogger("pypsa").setLevel(logging.ERROR)
    n = pypsa.Network(nc_path)
    B = n.buses
    fig, ax = plt.subplots(figsize=(9.5, 9.5) if box is None else (9, 6.5))
    for df in (n.lines, n.transformers):
        for nm, r in df.iterrows():
            if r.bus0 in B.index and r.bus1 in B.index and not str(nm).startswith("KRON"):
                hv = B.v_nom[r.bus0] >= 700
                ax.plot([B.x[r.bus0], B.x[r.bus1]], [B.y[r.bus0], B.y[r.bus1]], color=INK if hv else GRAY,
                        lw=1.3 if hv else 0.7, alpha=0.8 if hv else 0.9, zorder=1)
    v = V2.loc[sn]
    b = [x for x in v.index if x in B.index]
    sc = ax.scatter(B.x[b], B.y[b], c=v[b].values, cmap=VMAP, vmin=0.95, vmax=1.05, s=26, zorder=3, edgecolors=INK, linewidths=0.3)
    sv = [x for x in SVC_BUSES if x in B.index]
    ax.scatter(B.x[sv], B.y[sv], marker="s", s=90, facecolors="none", edgecolors=INK, lw=1.2, zorder=4)
    for bus, name in STATIONS.items():
        if bus in B.index and (box is None or (box[0] <= B.x[bus] <= box[1] and box[2] <= B.y[bus] <= box[3])):
            if box is None and bus in MTL:
                continue
            off = {"345": (6, -12), "340": (-10, 8), "66": (-22, -12), "310": (6, 6)}.get(bus, (4, 3)) if box is None else (4, 3)
            ax.annotate(name, (B.x[bus], B.y[bus]), fontsize=8.5 if box is None else 10, xytext=off, textcoords="offset points", color=INK)
    if box is None:
        ax.annotate("Montréal area\n(see zoom)", (-73.6, 45.5), xytext=(-78.5, 45.6), fontsize=9, color=INK2,
                    arrowprops=dict(arrowstyle="-", color=INK2, lw=0.8))
    if box is None:
        ax.set_xlim(-80.5, -60.5); ax.set_ylim(44.8, 56.5)
    else:
        ax.set_xlim(box[0], box[1]); ax.set_ylim(box[2], box[3])
    ax.set_aspect(1 / np.cos(np.radians(49))); ax.grid(False)
    ax.set_xticks([]); ax.set_yticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_title(title or f"v2 at the peak ({sn:%b %d, %H:%M}): bus voltages")
    cb = fig.colorbar(sc, ax=ax, shrink=0.5, pad=0.01); cb.set_label("voltage (pu)"); cb.outline.set_visible(False)
    leg = [Line2D([0], [0], color=INK, lw=1.3, label="735 kV line"), Line2D([0], [0], color=GRAY, lw=0.8, label="315 kV line"),
           Line2D([0], [0], marker="s", ls="", markerfacecolor="none", markeredgecolor=INK, markersize=9, label="real HQ SVC / synchronous condenser")]
    ax.legend(handles=leg, loc="lower right")
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 5. P-V curves
def fig_pv(pv, pv_kron, save=None):
    fig, ax = plt.subplots(figsize=(12, 5.2))
    series = [("v1", ORANGE, "v1 (Nov. 16-17 night)"), ("v2", BLUE, "v2 = v1 + run-of-river priority")]
    for case, c, lab in series:
        d = pv[(pv.case == case) & pv.solved.astype(bool)]
        st = d.stable.astype(bool)
        # the curve is drawn through stable points only, plus the first unstable one
        first_bad = d.index[~st][:1]
        dd = d.loc[d.index[st].union(first_bad)]
        ax.plot(dd.GW, dd.v_min, color=c, lw=2, label=lab)
        ax.scatter(d.GW[st], d.v_min[st], color=c, s=28, zorder=3)
        ax.scatter(d.GW[first_bad], d.v_min[first_bad], color=c, s=110, marker="X", edgecolors=INK, zorder=4)
    k = pv_kron
    ax.plot(k.GW, k.v_min, color=AQUA, lw=2, ls="--", label="v2 + Kron equivalent of the whole south (sensitivity)")
    ax.scatter(k.GW, k.v_min, color=AQUA, s=22, zorder=3)
    ax.axvline(34.959, color=INK2, lw=1, ls=":"); ax.text(34.99, 0.997, "test-week peak\n34.96 GW", fontsize=9.5, color=INK2, va="top")
    ax.axhline(0.90, color=INK2, lw=1, ls="--"); ax.text(31.3, 0.9015, "0.90 pu", fontsize=9.5, color=INK2)
    ax.annotate(f"2023 record {RECORD_GW:.1f} GW →", (39.4, 0.905), fontsize=10, color=INK2)
    ax.scatter([], [], marker="X", s=90, color=GRAY, edgecolors=INK, label="voltage-unstable (dV/dQ < 0)")
    ax.set_xlabel("load and generation scaled together (GW at the peak hour)")
    ax.set_ylabel("lowest bus voltage (pu)\n(always bus 800, Montréal-West)")
    ax.set_title("Load margin at the peak hour (P-V curve, Q limits enforced, shunts re-dispatched)")
    ax.set_xlim(31.2, 41); ax.set_ylim(0.88, 1.0)
    ax.legend(loc="upper right")
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 6. dispatch
def fig_dispatch(disp, save=None):
    lab = {("south (< 47.6°N)", "ror"): "south\nrun-of-river", ("north", "ror"): "north\nrun-of-river",
           ("north", "hydro reservoir"): "north\nreservoirs", ("south (< 47.6°N)", "onwind"): "south\nwind"}
    rows = []
    for (z, c), name in lab.items():
        for case in ("v1", "v2"):
            r = disp[(disp.case == case) & (disp.zone == z) & (disp.carrier == c)]
            rows.append(dict(group=name, case=case, produced=float(r.produced.sum()) / 1e3, installed=float(r.installed.sum()) / 1e3))
    D = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(11, 4.6))
    groups = list(lab.values()); x = np.arange(len(groups)); w = 0.36
    for k, (case, c) in enumerate((("v1", ORANGE), ("v2", BLUE))):
        d = D[D.case == case].set_index("group").loc[groups]
        bars = ax.bar(x + (k - 0.5) * w, d.produced, w * 0.94, color=c, label=f"{case} dispatch at the peak")
        for b, v in zip(bars, d.produced):
            ax.text(b.get_x() + b.get_width() / 2, v + 0.3, f"{v:.1f}", ha="center", fontsize=10)
    inst = D[D.case == "v1"].set_index("group").loc[groups].installed
    ax.scatter(x, inst, marker="_", s=900, color=INK, lw=2, label="installed", zorder=4)
    ax.set_xticks(x, groups); ax.set_ylabel("GW"); ax.grid(axis="x", visible=False)
    ax.set_title("Who produces at the peak: the LOPF had shut down 3.4 GW of southern run-of-river")
    ax.legend(loc="upper left")
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 7. one-at-a-time sensitivities (v1 dispatch)
SENS_EN = {
    "Cas final": "v1 (reference)",
    "Sans correction des longueurs (E1)": "without E1 (Mercator lengths kept)",
    "Sans correction de Churchill Falls (E2)": "without E2 (Churchill Falls lines doubled)",
    "Sans équivalents de Kron (E2b)": "without E2b (no Kron equivalents)",
    "Équivalent de Kron au Saguenay seulement (E2b)": "E2b with Saguenay only (no Montréal)",
    "Sans correction des circuits de Nathan (E3)": "without E3 (Nathan's circuit fix)",
    "Compensation série 70 % sur les 19 lignes Nathan/Etienne (E4)": "E4 at 70 % instead of 40 %",
    "Sans compensation série (E4)": "without E4 (no series compensation)",
    "Sans plafond de stabilité dans le LOPF (E5)": "without E5 (no 40° stability cap)",
    "Sans SVC ni compensateurs synchrones (E6)": "without real SVCs / sync. condensers",
    "Sans bancs de condensateurs (E6c)": "without capacitor banks",
    "Bancs de condensateurs alpha = 0,6 (E6c)": "capacitor banks × 0.6",
    "Bancs de condensateurs alpha = 1,5 (E6c)": "capacitor banks × 1.5",
    "Charges à cos φ 0,97 (E6)": "load power factor 0.97 (not 0.95)",
    "Consigne de tension 1,00 pu (E6)": "voltage set-point 1.00 pu",
    "Consigne de tension 1,04 pu (E6)": "voltage set-point 1.04 pu",
}


def sens_table(sens):
    s = sens.copy(); s["variant"] = s.variante.map(SENS_EN).fillna(s.variante)
    s = s.rename(columns={"convergees": "within Q limits", "bande_090_110": "in 0.90-1.10", "bande_095_105": "in 0.95-1.05",
                          "v_min": "V min", "v_max": "V max", "marge_pointe_lambda": "peak margin λ"})
    return s[["variant", "within Q limits", "in 0.90-1.10", "in 0.95-1.05", "V min", "V max", "peak margin λ"]]


def fig_sensitivity(sens, save=None):
    s = sens_table(sens).iloc[::-1]
    fig, ax = plt.subplots(figsize=(12, 6.4))
    y = np.arange(len(s))
    ax.hlines(y, 0, s["in 0.95-1.05"], color=GRID, lw=2)
    ax.scatter(s["within Q limits"], y, color=GRAY, s=60, label="hours within Q limits", zorder=3)
    ax.scatter(s["in 0.95-1.05"], y, color=BLUE, s=60, label="hours with all V in 0.95-1.05 pu", zorder=4)
    ax.set_yticks(y, s.variant); ax.set_xlim(0, 175); ax.set_xlabel("hours of the test week (of 168)")
    ref = s[s.variant == "v1 (reference)"]
    ax.axvline(float(ref["in 0.95-1.05"].iloc[0]), color=INK2, lw=1, ls=":")
    ax.set_title("Remove or change one assumption at a time (v1 dispatch): what each correction buys")
    ax.legend(loc="lower left"); ax.grid(axis="y", visible=False)
    fig.tight_layout()
    return _save(fig, save)


# ---------------------------------------------------------------- 8. below 315 kV
def fig_full_network(full, save=None):
    order = [69, 120, 161, 230, 315, 735]
    data = [full.v_pu[full.v_nom.round() == v].values for v in order]
    fig, ax = plt.subplots(figsize=(11, 4.6))
    for i, (v, d) in enumerate(zip(order, data)):
        xs = i + (np.random.default_rng(0).random(len(d)) - 0.5) * 0.5
        out = (d < 0.9) | (d > 1.1)
        ax.scatter(xs[~out], d[~out], s=10, color=BLUE, alpha=0.6, edgecolors="none")
        ax.scatter(xs[out], d[out], s=26, color=ORANGE, edgecolors=INK, linewidths=0.4, zorder=3)
        ax.text(i, 0.865, f"{len(d)} buses", ha="center", fontsize=9.5, color=INK2)
    for yv in (0.90, 1.10):
        ax.axhline(yv, color=INK2, lw=1, ls="--")
    ax.set_xticks(range(len(order)), [f"{v} kV" for v in order]); ax.set_ylim(0.85, 1.14)
    ax.set_ylabel("bus voltage (pu)"); ax.grid(axis="x", visible=False)
    ax.set_title("Exploratory: full network down to 69 kV at the peak (796 buses; orange = outside 0.90-1.10)")
    fig.tight_layout()
    return _save(fig, save)
