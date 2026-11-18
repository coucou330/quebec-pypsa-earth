"""Réseau complet du Québec (toutes tensions) pour l'AC PF — exploratoire.

Principe : on part du modèle 315 kV construit par reseau315.build() SANS équivalents de Kron (mêmes
corrections E1-E6, même LOPF), et on remet dessous le vrai réseau < 315 kV d'elec_full.nc :
  - bus < 315 kV du Québec = bus d'origine des charges du modèle (la réduction de Nathan avait déplacé
    539 charges sur 647 vers le bus >= 315 kV le plus proche) + bus voisins à moins de `rayon_km` ;
  - lignes < 315 kV : paramètres des types PyPSA-Earth, corrigés comme dans reseau315 (longueur : règle
    hybride ; types 50 Hz : x et b × 1,2) ;
  - transformateurs d'elec_full (génériques : 2 000 MVA, x = 0,1 pu, r mis à x / 50) ;
  - charges remises sur leur bus d'origine ; bancs de condensateurs (alpha × Q de la charge à la pointe)
    replacés aux nouveaux bus de charge ;
  - centrales, régulateurs, SVC/CS, inductances : laissés où ils sont dans le modèle 315 kV.
Rien n'est écrit dans network/networks_current/."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
from _commun import *
import scipy.sparse as sp, scipy.sparse.linalg as spl
from scipy.spatial import cKDTree


class SparseEngine(R.ACEngine):
    """Même moteur que reseau315.ACEngine, avec un Newton creux (le réseau complet a ~5 fois plus de bus)."""
    def newton(self, Y, S, V0, pv, pq, tol=1e-3, maxit=30):
        Ys = sp.csr_matrix(Y); V = V0.copy(); pvpq = np.r_[pv, pq].astype(int); npv = len(pvpq); pq = np.asarray(pq, int)
        for it in range(maxit + 1):
            I = Ys @ V; mis = V * np.conj(I) - S; F = np.r_[mis[pvpq].real, mis[pq].imag]; err = np.abs(F).max()
            if not np.isfinite(err): return V, False, it
            if err < tol: return V, True, it
            dV = sp.diags(V); dI = sp.diags(I); dVn = sp.diags(V / np.abs(V))
            dth = (1j * dV @ (dI - Ys @ dV).conj()).tocsr(); dvm = (dV @ (Ys @ dVn).conj() + dI.conj() @ dVn).tocsr()
            J = sp.bmat([[dth[pvpq][:, pvpq].real, dvm[pvpq][:, pq].real], [dth[pq][:, pvpq].imag, dvm[pq][:, pq].imag]], format="csc")
            try: dx = spl.spsolve(J, -F)
            except Exception: return V, False, it
            if not np.isfinite(dx).all(): return V, False, it
            th0, vm0 = np.angle(V), np.abs(V)
            for a in (1.0, 0.5, 0.25, 0.1):
                th = th0.copy(); vm = vm0.copy(); th[pvpq] += a * dx[:npv]; vm[pq] += a * dx[npv:]
                Vt = vm * np.exp(1j * th); mt = Vt * np.conj(Ys @ Vt) - S; Ft = np.r_[mt[pvpq].real, mt[pq].imag]
                if np.isfinite(Ft).all() and np.abs(Ft).max() < err or a == 0.1: V = Vt; break
            if np.abs(V).min() < 0.3 or np.abs(V).max() > 2.0: return V, False, it
        return V, False, maxit


def build_complet(full, T, rayon_km=25.0, alpha=1.0, beta=R.BETA_50HZ, build_kw=None, lopf=None, log=None, max_ratio=0.10, exclure=()):
    log = [] if log is None else log
    lopf0 = R.lopf
    if lopf is not None: R.lopf = lopf
    try:
        with R.quiet():
            n, info = R.build(T=T, kron=False, alpha_banks=alpha, **(build_kw or {}))
    finally:
        R.lopf = lopf0
    backbone = set(n.buses.index)
    FB = full.buses
    orig = full.loads.bus.reindex(n.loads.index)                       # bus d'origine de chaque charge
    seed = set(orig.dropna().values) | (backbone & set(FB.index))
    # candidats : bus AC < 315 kV d'elec_full à moins de rayon_km d'un bus de charge du Québec ou d'un bus du modèle
    cand = FB[(FB.carrier == "AC") & (FB.v_nom < 315) & ~FB.index.isin(backbone)]
    sx = np.cos(np.radians(48.0)) * 111.2; sy = 111.2
    S = FB.loc[[b for b in seed if b in FB.index]]
    tree = cKDTree(np.c_[S.x * sx, S.y * sy]); d, _ = tree.query(np.c_[cand.x * sx, cand.y * sy])
    # on ne garde un bus sans charge du Québec que s'il est près d'un bus connu ; les bus de charge sont toujours gardés
    local = set(cand.index[d <= rayon_km]) | (set(orig.dropna().values) - backbone)
    local = [b for b in local if b in FB.index and FB.at[b, "carrier"] == "AC"]
    with R.quiet():
        full.calculate_dependent_values()
        n.add("Bus", local, v_nom=FB.v_nom[local].values, x=FB.x[local].values, y=FB.y[local].values, carrier="AC", v_mag_pu_set=1.0)
    keep = backbone | set(local)
    sm = float(T.loc[T.reliable, "sinuosity"].median())
    L = full.lines[full.lines.bus0.isin(keep) & full.lines.bus1.isin(keep)]
    L = L[~(L.bus0.isin(backbone) & L.bus1.isin(backbone))]
    f = pd.Series(1.0, index=L.index); inT = L.index.intersection(T.index)
    f[inT] = np.where(T.loc[inT, "reliable"], T.loc[inT, "true_km"], sm * T.loc[inT, "gc_km"]) / T.loc[inT, "length_stored"]
    fb = np.where(L.v_nom < 315, beta, 1.0)
    names = ["BT " + l for l in L.index]
    with R.quiet():
        n.add("Line", names, bus0=L.bus0.values, bus1=L.bus1.values, r=(L.r * f).values, x=(L.x * f * fb).values, b=(L.b * f * fb).values,
              length=(L.length * f).values, num_parallel=1.0, s_nom=5000.0, v_nom=L.v_nom.values, carrier="AC")
    Tr = full.transformers[full.transformers.bus0.isin(keep) & full.transformers.bus1.isin(keep)]
    Tr = Tr[~(Tr.bus0.isin(backbone) & Tr.bus1.isin(backbone))]
    with R.quiet():
        n.add("Transformer", ["BT " + t for t in Tr.index], bus0=Tr.bus0.values, bus1=Tr.bus1.values, x=Tr.x.values, r=(Tr.x / 50.0).values,
              s_nom=Tr.s_nom.values, model="t")
        n.calculate_dependent_values(); n.determine_network_topology()
    slack_bus = n.generators.bus[n.generators.control == "Slack"]
    main = n.buses.sub_network[[b for b in slack_bus if n.buses.at[b, "carrier"] == "AC"]].value_counts().idxmax()
    inmain = set(n.buses.index[n.buses.sub_network == main])
    # charges : retour au bus d'origine s'il est relié au réseau principal
    back = [l for l in n.loads.index if isinstance(orig.get(l), str) and orig[l] in inmain and orig[l] != n.loads.at[l, "bus"]]
    stay = [l for l in n.loads.index if isinstance(orig.get(l), str) and orig[l] not in backbone and orig[l] not in inmain]
    # charges que le réseau local ne peut physiquement pas porter : elles restent au bus >= 315 kV.
    # Critère : charge maximale du bus > max_ratio × somme des V^2/x de ses branches (ex. bus 3559, 69 kV :
    # 851 MW au bout d'une seule ligne 69 kV dont V^2/x = 835 MW ; c'est un artefact de la répartition régionale
    # de la demande, rescale_demand_regional.py, pas une vraie charge).
    force = pd.Series(0.0, index=n.buses.index)
    for df, z in ((n.lines, n.lines.x / n.lines.v_nom ** 2), (n.transformers, n.transformers.x / n.transformers.s_nom)):
        y = 1.0 / z.abs().replace(0, np.nan)
        force = force.add(y.groupby(df.bus0).sum(), fill_value=0).add(y.groupby(df.bus1).sum(), fill_value=0)
    lmax = n.loads_t.p_set[back].max().groupby(orig[back].values).sum()
    trop = set(lmax.index[lmax > max_ratio * force.reindex(lmax.index).fillna(0)]) | (set(exclure) & set(lmax.index))   # exclure : bus dont les charges restent au 315 kV (c5_elagage.py)
    trop_l = [l for l in back if orig[l] in trop]
    back = [l for l in back if orig[l] not in trop]
    stay = stay + trop_l
    p_back = n.loads_t.p_set[back].loc[PEAK].sum(); p_stay = n.loads_t.p_set[stay].loc[PEAK].sum()
    n.loads.loc[back, "bus"] = orig[back].values
    # bus < 315 kV non reliés au réseau principal : retirés (avec leurs branches)
    drop = [b for b in local if b not in inmain]
    with R.quiet():
        n.remove("Line", n.lines.index[n.lines.bus0.isin(drop) | n.lines.bus1.isin(drop)])
        n.remove("Transformer", n.transformers.index[n.transformers.bus0.isin(drop) | n.transformers.bus1.isin(drop)])
        n.remove("Bus", drop)
        # condensateurs : replacés aux nouveaux bus de charge
        n.remove("ShuntImpedance", n.shunt_impedances.index[n.shunt_impedances.kind == "capacitor"])
    if alpha:
        banks = R.capacitor_banks_from_load(n, alpha=alpha, min_q=0.5)
        R.add_capacitor_banks(n, banks)
    with R.quiet():
        n.calculate_dependent_values(); n.determine_network_topology()
    info["complet"] = dict(bus_ajoutes=len(local) - len(drop), bus_non_relies=len(drop), lignes=len(names), transfos=len(Tr),
                           charges_redescendues=len(back), MW_redescendus_pointe=float(p_back), charges_restees_au_315=len(stay), MW_restes=float(p_stay), bus_surcharges=sorted(trop),
                           lignes_sans_geometrie=int(len(L) - len(inT)), km_lignes_bt=float((L.length * f).sum()))
    return n, info


def bilan(n, eng, r):
    vm = pd.Series(np.abs(r["V"]), index=eng.buses); vn = n.buses.v_nom.reindex(eng.buses)
    cls = pd.cut(vn, [0, 100, 200, 300, 400, 800], labels=["< 100 kV", "100-200 kV", "200-300 kV", "315-345 kV", "735 kV"])
    return vm.groupby(cls, observed=True).agg(["size", "min", "median", "max"]).round(3), vm
