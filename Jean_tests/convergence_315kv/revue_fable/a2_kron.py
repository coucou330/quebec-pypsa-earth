"""A2 : vérification numérique de la réduction de Kron (boîtes de Montréal et du Saguenay).
(a) AC PF du sous-réseau complet (bus frontières + bus éliminés, injections nulles aux bus éliminés) contre
    AC PF du réseau réduit (Yeq) : mêmes tensions aux bus frontières ?
(b) ce qui est réellement ajouté au modèle PyPSA (add_kron_equivalent : branches |y| < 1e-5 jetées, r borné
    à 1e-6, partie réelle des shunts jetée) reproduit-il Yeq ?
(c) hypothèses de la fonction : branches frontière-frontière « déjà dans le modèle » ? règle de longueur ?"""
from _commun import *
import inspect
full = R.load(os.path.join(R.CURRENT, "elec_full.nc"))
T = R.mercator_table(full)
n, eng = final()
src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))

def full_Y(box):
    """Recopie de la construction de Y dans kron_equivalent (mêmes paramètres), sans la réduction."""
    B = full.buses
    region = set(B.index[(B.y.between(*box["lat"])) & (B.x.between(*box["lon"]))])
    mb = set(src.buses.index)
    elim = sorted(b for b in region if b not in mb); bound = sorted(b for b in region if b in mb)
    idx = {b: i for i, b in enumerate(bound + elim)}; N = len(idx); Y = np.zeros((N, N), complex)
    full.calculate_dependent_values()
    skipped = []
    for name, r in full.lines[full.lines.bus0.isin(region) & full.lines.bus1.isin(region)].iterrows():
        a, b = r.bus0, r.bus1
        if a in mb and b in mb:
            skipped.append(("Line", name, a, b, r.v_nom)); continue
        f = T.at[name, "factor"] if name in T.index else 1.0
        fb = R.BETA_50HZ if r.v_nom < 315 else 1.0
        z = (r.r * f + 1j * r.x * f * fb) / r.v_nom ** 2
        if abs(z) == 0: continue
        ysh = 1j * r.b * f * fb * r.v_nom ** 2 / 2
        i, j = idx[a], idx[b]
        Y[i, i] += 1 / z + ysh; Y[j, j] += 1 / z + ysh; Y[i, j] -= 1 / z; Y[j, i] -= 1 / z
    for name, r in full.transformers[full.transformers.bus0.isin(region) & full.transformers.bus1.isin(region)].iterrows():
        a, b = r.bus0, r.bus1
        if a in mb and b in mb:
            skipped.append(("Transformer", name, a, b, 0)); continue
        z = (r.r + 1j * r.x) / r.s_nom
        if abs(z) == 0: z = 1j * 0.1 / r.s_nom
        i, j = idx[a], idx[b]
        Y[i, i] += 1 / z; Y[j, j] += 1 / z; Y[i, j] -= 1 / z; Y[j, i] -= 1 / z
    return Y, bound, elim, skipped

def nr(Y, S, slack, V0=None, tol=1e-9):
    N = len(S); pq = np.array([i for i in range(N) if i != slack]); V = np.ones(N, complex) if V0 is None else V0.copy()
    for it in range(50):
        mis = V * np.conj(Y @ V) - S; F = np.r_[mis[pq].real, mis[pq].imag]
        if np.abs(F).max() < tol: return V, True, it
        I = Y @ V; dV = np.diag(V); dVn = np.diag(V / np.abs(V))
        dth = 1j * dV @ np.conj(np.diag(I) - Y @ dV); dvm = dV @ np.conj(Y @ dVn) + np.conj(np.diag(I)) @ dVn
        J = np.block([[dth[np.ix_(pq, pq)].real, dvm[np.ix_(pq, pq)].real], [dth[np.ix_(pq, pq)].imag, dvm[np.ix_(pq, pq)].imag]])
        dx = np.linalg.solve(J, -F); th = np.angle(V); vm = np.abs(V)
        th[pq] += dx[:len(pq)]; vm[pq] += dx[len(pq):]; V = vm * np.exp(1j * th)
    return V, False, it

rng = np.random.default_rng(0)
for reg, box in R.KRON_REGIONS.items():
    print(f"\n===== {reg}")
    Y, bound, elim, skipped = full_Y(box)
    nbd = len(bound)
    keep = np.r_[np.arange(nbd), nbd + np.where(np.abs(Y[nbd:, nbd:]).sum(axis=1) > 0)[0]]
    Yk = Y[np.ix_(keep, keep)]
    Yeq = Yk[:nbd, :nbd] - Yk[:nbd, nbd:] @ np.linalg.solve(Yk[nbd:, nbd:], Yk[nbd:, :nbd])
    br, sh, info = R.kron_equivalent(full, src, T, box=box)
    print(f"bus frontières {nbd}, éliminés {len(elim)} (dont {len(elim) - (len(keep) - nbd)} isolés), branches équivalentes {len(br)}")
    # composantes connexes du sous-réseau : un bus frontière sans lien avec les autres ne peut pas être dans une PF
    import scipy.sparse.csgraph as cg
    nc, lab = cg.connected_components(np.abs(Yk) > 0)
    comp_b = pd.Series(lab[:nbd], index=bound)
    print(f"composantes connexes du sous-réseau régional (sans les branches frontière-frontière) : {nc} ; frontières par composante : {comp_b.value_counts().to_dict()}")
    big = comp_b.value_counts().idxmax(); sel = np.where(lab == big)[0]; selb = [i for i in sel if i < nbd]
    # (a) AC PF complet contre réduit, sur la plus grande composante
    Ya = Yk[np.ix_(sel, sel)]; nb_ = len(selb)
    Yr = Yeq[np.ix_(selb, selb)]
    S = np.zeros(len(sel), complex)
    S[1:nb_] = rng.uniform(-150, 150, nb_ - 1) + 1j * rng.uniform(-50, 50, nb_ - 1)   # MW / MVAr aux bus frontières (slack = le 1er)
    Vf, okf, itf = nr(Ya, S, 0); Vr, okr, itr = nr(Yr, S[:nb_], 0)
    print(f"(a) AC PF complet ({len(sel)} bus) : converge {okf} en {itf} it ; réduit ({nb_} bus) : {okr} en {itr} it")
    print(f"    écart max de tension aux bus frontières : |V| {np.abs(np.abs(Vf[:nb_]) - np.abs(Vr)).max():.2e} pu, angle {np.degrees(np.abs(np.angle(Vf[:nb_]) - np.angle(Vr))).max():.2e}°"
          f"   (|V| entre {np.abs(Vr).min():.3f} et {np.abs(Vr).max():.3f})")
    # (b) ce que PyPSA a réellement dans le réseau final
    ix = {b: i for i, b in enumerate(bound)}; Yimp = np.zeros((nbd, nbd), complex)
    pre = f"KRON {reg}"
    for nm, r in n.lines[n.lines.index.str.startswith(pre)].iterrows():
        y = 1 / ((r.r + 1j * r.x) / r.v_nom ** 2); i, j = ix[r.bus0], ix[r.bus1]
        Yimp[i, i] += y; Yimp[j, j] += y; Yimp[i, j] -= y; Yimp[j, i] -= y
    for nm, r in n.transformers[n.transformers.index.str.startswith(pre)].iterrows():
        y = 1 / ((r.r + 1j * r.x) / r.s_nom); i, j = ix[r.bus0], ix[r.bus1]
        Yimp[i, i] += y; Yimp[j, j] += y; Yimp[i, j] -= y; Yimp[j, i] -= y
    for nm, r in n.shunt_impedances[n.shunt_impedances.index.str.startswith(pre + " sh")].iterrows():
        Yimp[ix[r.bus], ix[r.bus]] += 1j * r.b * n.buses.v_nom[r.bus] ** 2
    D = Yimp - Yeq
    off = ~np.eye(nbd, dtype=bool)
    print(f"(b) implémenté contre Yeq : écart hors diagonale max {np.abs(D[off]).max():.2e} pu (|Yeq| hors diag. max {np.abs(Yeq[off]).max():.2e}) ; "
          f"somme des |y| jetés (< 1e-5) {np.abs(Yeq[off])[np.abs(Yeq[off]) <= 1e-5].sum() / 2:.2e}")
    rowsum = Yeq.sum(axis=1)
    print(f"    shunts équivalents : partie réelle jetée, total {rowsum.real.sum():.4f} pu (= {rowsum.real.sum():.1f} MW à 1 pu) ; partie imaginaire gardée {rowsum.imag.sum():.1f} MVAr")
    z = br.z_pu.values
    print(f"    branches à r < 0 (borné à 1e-6) : {(z.real < 0).sum()}/{len(z)} ; à x < 0 : {(z.imag < 0).sum()}/{len(z)} ; x min {z.imag.min():.2e}, x max {z.imag.max():.2e} pu")
    Vr2, ok2, _ = nr(Yimp[np.ix_(selb, selb)], S[:nb_], 0)
    print(f"    AC PF avec l'équivalent tel qu'implémenté : converge {ok2} ; écart |V| max aux frontières contre le réseau complet {np.abs(np.abs(Vf[:nb_]) - np.abs(Vr2)).max():.2e} pu, angle {np.degrees(np.abs(np.angle(Vf[:nb_]) - np.angle(Vr2))).max():.2e}°")
    # (c) branches frontière-frontière sautées : sont-elles dans le modèle ?
    miss = []
    for kind, name, a, b, v in skipped:
        df = src.lines if kind == "Line" else src.transformers
        if name not in df.index: miss.append((kind, name, a, b, v))
    print(f"(c) branches frontière-frontière sautées : {len(skipped)} ; absentes du modèle (donc perdues) : {len(miss)} {miss[:8]}")
    # règle de longueur utilisée dans Kron : T.factor (géodésique/Mercator), pas la règle hybride
    Lr = full.lines[full.lines.bus0.isin(bound + elim) & full.lines.bus1.isin(bound + elim)]
    Lr = Lr[~(Lr.bus0.isin(bound) & Lr.bus1.isin(bound))]
    inT = Lr.index.intersection(T.index)
    hyb = np.where(T.loc[inT, "reliable"], T.loc[inT, "true_km"], 1.084 * T.loc[inT, "gc_km"]) / T.loc[inT, "length_stored"]
    used = T.loc[inT, "factor"].values
    w = Lr.loc[inT, "x"].values
    print(f"    longueurs : {len(Lr)} lignes éliminées, {len(Lr) - len(inT)} sans géométrie (facteur 1, donc encore en Mercator ?) ; facteur utilisé (géod./Mercator) médian {np.median(used):.3f} "
          f"contre règle hybride {np.median(hyb):.3f} ; écart pondéré par x : {((hyb - used) * w).sum() / (used * w).sum():+.1%} ; tracés non fiables {int((~T.loc[inT, 'reliable']).sum())}")
