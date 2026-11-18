"""Outils de la partie B (aucune modification de reseau315.py).
ATTENTION : continuation(), opf_*() et newton_statisme()/marge_statisme() sont des essais ABANDONNÉS
(cf. RAPPORT.md, B1) : la continuation trouve de fausses solutions (régulateurs bloqués à Qmax à 1,2 pu),
l'OPF de pandapower ne converge sur aucun cas, le Newton à statisme ne converge pas près du nez.
Seuls Ymat() et slack_reparti() sont utiles ; le critère retenu est celui de b1d_stabilite_vq.py.
Description d'origine :
- solve_from / continuation : Newton avec limites Q, en repartant de la solution ET de l'état PV/PQ du pas
  précédent (une vraie continuation en charge, par petits pas), au lieu de repartir à froid à chaque pas ;
- opf_margin : existence d'un point de fonctionnement admissible (OPF AC de pandapower : limites Q, bande de
  tension, consignes de tension libres), indépendante de la logique PV/PQ ;
- slack_reparti : pertes réparties sur les centrales au lieu d'être toutes produites par LG-2."""
from _commun import *


def solve_from(eng, Y, P, Qf, V, state, max_rounds=80, max_flips=4, one_at_a_time=True):
    """Comme ACEngine.solve, mais part d'un état PV/PQ donné et ne bascule qu'un régulateur par tour."""
    nreg = len(eng.regs); state = state.copy(); flips = np.zeros(nreg, int); V_ok = None; st_ok = state.copy()
    for rnd in range(max_rounds):
        pv = sorted({eng.reg_bus[j] for j in range(nreg) if state[j] == 0})
        pvs = set(pv); pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]
        Qs = Qf.copy()
        for j in np.where(state != 0)[0]:
            Qs[eng.reg_bus[j]] += eng.qmax[j] if state[j] > 0 else eng.qmin[j]
        Vm = np.abs(V); Vm[pv] = eng.vset[pv]; Vm[eng.sl] = eng.vset[eng.sl]
        V1, ok, it = eng.newton(Y, P + 1j * Qs, Vm * np.exp(1j * np.angle(V)), np.array(pv, int), np.array(pq, int))
        if not ok:
            break
        V = V1; V_ok, st_ok = V.copy(), state.copy()
        q = (V * np.conj(Y @ V)).imag - Qf
        viol, back = [], []
        for j, b in enumerate(eng.reg_bus):
            if state[j] == 0:
                if q[b] > eng.qmax[j] + 1: viol.append((q[b] - eng.qmax[j], j, 1))
                elif q[b] < eng.qmin[j] - 1: viol.append((eng.qmin[j] - q[b], j, -1))
            elif flips[j] < max_flips:
                d = abs(V[b]) - eng.vset[b]
                if (state[j] > 0 and d > 1e-4) or (state[j] < 0 and d < -1e-4): back.append((abs(d), j))
        if viol:
            for _, j, s in sorted(viol, reverse=True)[:1 if one_at_a_time else 3]:
                state[j] = s; flips[j] += 1
        elif back:
            _, j = max(back); state[j] = 0; flips[j] += 1
        else:
            break
    if V_ok is None:
        return dict(ok=False)
    V = V_ok; q = (V * np.conj(Y @ V)).imag - Qf
    nv = sum(1 for j, b in enumerate(eng.reg_bus) if st_ok[j] == 0 and (q[b] > eng.qmax[j] + 1 or q[b] < eng.qmin[j] - 1))
    comp = max([abs(V[b]) - eng.vset[b] for j, b in enumerate(eng.reg_bus) if st_ok[j] > 0] + [0.0])
    return dict(ok=True, V=V, state=st_ok, n_viol=nv, compl=comp, p_slack=(V * np.conj(Y @ V)).real[eng.sl] - P[eng.sl])


def Ymat(eng, k, c):
    return eng.Y0 + np.diag(1j * (np.broadcast_to(k, (eng.nb,)) * eng.b_react + np.broadcast_to(c, (eng.nb,)) * eng.b_cap))


def continuation(eng, sn, k, c, V0=None, state0=None, step=0.0025, lam_max=1.4, P0=None, Q0=None, verbose=False):
    """Marge de charge par continuation : charge et production × lambda, shunts figés (k, c), départ de chaque
    pas = solution et état du pas précédent. S'arrête au premier pas sans solution dans les limites Q.
    Retourne un DataFrame (lam, v_min, bus, n_pq, compl) ; compl = pire dépassement de consigne d'un régulateur en butée haute."""
    P0 = eng.P.loc[sn].values.copy() if P0 is None else P0; Q0 = eng.Q.loc[sn].values.copy() if Q0 is None else Q0
    Y = Ymat(eng, k, c)
    if V0 is None:
        r = eng.solve(sn, k=k, c=c); V0, state0 = r["V"], r["state"]
    V, st = V0.copy(), state0.copy(); rows = []
    for lam in np.round(np.arange(1.0, lam_max + 1e-9, step), 5):
        r = solve_from(eng, Y, P0 * lam, Q0 * lam, V, st)
        if not r["ok"] or r["n_viol"] > 0:
            # second essai : repartir de tous les régulateurs en PV, depuis la dernière solution
            r = solve_from(eng, Y, P0 * lam, Q0 * lam, V, np.zeros(len(st), int))
            if not r["ok"] or r["n_viol"] > 0:
                break
        V, st = r["V"], r["state"]; vm = np.abs(V)
        rows.append(dict(lam=lam, v_min=vm.min(), bus=eng.buses[vm.argmin()], n_pq=int((st != 0).sum()), compl=r["compl"], p_slack=r["p_slack"]))
        if verbose: print(rows[-1], flush=True)
    return pd.DataFrame(rows)


def slack_reparti(eng, sn, k, c, frac=None, iters=6):
    """Répartit les pertes AC sur les centrales (au prorata de leur production, comme un réglage
    fréquence-puissance), au lieu de les faire toutes produire par le slack (LG-2). Retourne (P modifié, r).
    Le slack ne garde que sa part : on vise p_slack = dispatch LOPF du slack (0 ici) + sa part des pertes."""
    n = eng.n; P0 = eng.P.loc[sn].values.copy(); Qf = eng.Q.loc[sn].values
    gen = np.zeros(eng.nb)   # production par bus (hors slack) : centrales pilotables seulement
    su = n.storage_units; sp = n.storage_units_t.p_set.loc[sn].reindex(su.index).fillna(0.0)
    for nm, p in sp.items():
        b = su.at[nm, "bus"]
        if b in eng.ix and p > 0: gen[eng.ix[b]] += p
    g = n.generators; gp = n.generators_t.p_set.loc[sn].reindex(g.index).fillna(0.0)
    for nm, p in gp.items():
        if g.at[nm, "carrier"] in ("hydro", "ror", "OCGT") and g.at[nm, "bus"] in eng.ix and p > 0: gen[eng.ix[g.at[nm, "bus"]]] += p
    part = gen / gen.sum()           # participation (le bus slack compte avec ses centrales)
    Y = Ymat(eng, k, c); P = P0.copy(); r = None
    for it in range(iters):
        saveP = eng.P.loc[sn].copy(); eng.P.loc[sn] = P
        r = eng.solve(sn, k=k, c=c); eng.P.loc[sn] = saveP
        if not r["ok"]: return P, r
        extra = r["p_slack"]         # ce que le slack ajoute encore
        if abs(extra * (1 - part[eng.sl])) < 1.0: break
        add = extra * part; add[eng.sl] = 0.0
        P = P + add
    return P, r


def opf_net(n, eng, sn, k, c, lam=1.0, vmin=0.90, vmax=1.10, reactors_free=False, vset_free=True, P=None, Q=None):
    """Réseau pandapower pour une OPF AC de faisabilité : régulateurs = générateurs P = 0 avec limites Q,
    tensions libres dans [vmin, vmax] ; charges et centrales fixes (× lam) ; slack libre ; on minimise P slack.
    reactors_free : les inductances de ligne deviennent des sources de Q continues dans [-Q(k=1), 0]."""
    import pandapower as pp
    P = (eng.P.loc[sn].values if P is None else P) * lam; Q = (eng.Q.loc[sn].values if Q is None else Q) * lam
    net = pp.create_empty_network(sn_mva=100.0, f_hz=60)
    vn = eng.vnom
    for i, b in enumerate(eng.buses):
        pp.create_bus(net, vn_kv=vn[i], name=b, min_vm_pu=vmin, max_vm_pu=vmax)
    # branches : on passe par des impédances pu (base 100 MVA) -> identique à la matrice Y du moteur
    L = n.lines
    for nm, r in L.iterrows():
        if r.bus0 in eng.ix and r.bus1 in eng.ix:
            zb = r.v_nom ** 2 / 100.0
            i, j = eng.ix[r.bus0], eng.ix[r.bus1]
            pp.create_impedance(net, i, j, rft_pu=r.r / zb, xft_pu=r.x / zb, sn_mva=100.0)
            if r.b > 0:
                for e in (i, j): pp.create_shunt(net, e, q_mvar=-r.b / 2 * r.v_nom ** 2, p_mw=0.0, vn_kv=vn[e])
    for nm, r in n.transformers.iterrows():
        if r.bus0 in eng.ix and r.bus1 in eng.ix:
            sc = 100.0 / r.s_nom
            pp.create_impedance(net, eng.ix[r.bus0], eng.ix[r.bus1], rft_pu=r.r * sc, xft_pu=r.x * sc, sn_mva=100.0)
    sh = n.shunt_impedances
    for s, r in sh.iterrows():
        if r.bus in eng.ix and s not in eng.sw and r.b != 0:
            pp.create_shunt(net, eng.ix[r.bus], q_mvar=-r.b * n.buses.v_nom[r.bus] ** 2, p_mw=0.0, vn_kv=n.buses.v_nom[r.bus])
    kv = np.broadcast_to(k, (eng.nb,)); cv = np.broadcast_to(c, (eng.nb,))
    for i in range(eng.nb):
        bq = cv[i] * eng.b_cap[i] + (0.0 if reactors_free else kv[i] * eng.b_react[i])
        if bq != 0: pp.create_shunt(net, i, q_mvar=-bq, p_mw=0.0, vn_kv=vn[i])
        if reactors_free and eng.b_react[i] != 0:
            pp.create_sgen(net, i, p_mw=0.0, q_mvar=0.0, controllable=True, min_p_mw=0.0, max_p_mw=0.0,
                           min_q_mvar=eng.b_react[i] * 1.0, max_q_mvar=0.0, name="REACT")
    for i in range(eng.nb):
        if i != eng.sl and (P[i] != 0 or Q[i] != 0):
            pp.create_load(net, i, p_mw=-P[i], q_mvar=-Q[i], controllable=False)
    for j, b in enumerate(eng.reg_bus):
        lo, hi = (vmin, vmax) if vset_free else (eng.vset[b], eng.vset[b])
        pp.create_gen(net, b, p_mw=0.0, vm_pu=eng.vset[b], min_q_mvar=eng.qmin[j], max_q_mvar=eng.qmax[j], min_p_mw=0.0, max_p_mw=0.0, controllable=True, name=eng.regs[j])
        if not vset_free:
            net.bus.loc[b, ["min_vm_pu", "max_vm_pu"]] = [eng.vset[b] - 1e-4, eng.vset[b] + 1e-4]
    g = n.generators; sl = g.index[(g.control == "Slack") & (g.bus == eng.slack_bus)]
    qlo, qhi = float(g.loc[sl, "q_min"].sum()), float(g.loc[sl, "q_max"].sum())
    pp.create_ext_grid(net, eng.sl, vm_pu=eng.vset[eng.sl], min_p_mw=-1e5, max_p_mw=1e5, min_q_mvar=qlo, max_q_mvar=qhi)
    if P[eng.sl] != 0 or Q[eng.sl] != 0:
        pp.create_load(net, eng.sl, p_mw=-P[eng.sl], q_mvar=-Q[eng.sl], controllable=False)
    pp.create_poly_cost(net, 0, "ext_grid", cp1_eur_per_mw=1.0)
    return net


def opf_feasible(n, eng, sn, k, c, lam, **kw):
    import pandapower as pp
    net = opf_net(n, eng, sn, k, c, lam=lam, **kw)
    try:
        pp.runopp(net, init="pf" if False else "flat", calculate_voltage_angles=True, numba=False, PDIPM_MAX_IT=300)
        ok = bool(net.OPF_converged)
    except Exception as e:
        return dict(ok=False, err=type(e).__name__)
    if not ok: return dict(ok=False)
    return dict(ok=True, v_min=float(net.res_bus.vm_pu.min()), v_max=float(net.res_bus.vm_pu.max()), p_slack=float(net.res_ext_grid.p_mw.sum()),
                q_gen=net.res_gen.q_mvar.values.copy(), vm=net.res_bus.vm_pu.values.copy(), net=net)


def opf_margin(n, eng, sn, k, c, lams, **kw):
    rows = []
    for lam in lams:
        r = opf_feasible(n, eng, sn, k, c, lam, **kw)
        rows.append(dict(lam=lam, ok=r["ok"], v_min=r.get("v_min"), v_max=r.get("v_max"), p_slack=r.get("p_slack")))
        if not r["ok"]: break
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------------------------------------
# Régulateurs à statisme (caractéristique réelle d'un SVC / d'un régulateur de tension) :
#   Q = borne( (consigne - V) / statisme × Qmax , Qmin, Qmax )
# La complémentarité est alors automatique (Q = Qmax seulement si V est sous la consigne) et il n'y a plus
# de bascule PV/PQ à gérer : un seul Newton (semi-lisse). Statisme -> 0 redonne le PV/PQ idéal.
def newton_statisme(eng, Y, P, Qf, V0, slope=0.02, tol=1e-3, maxit=200):
    ns = np.array([i for i in range(eng.nb) if i != eng.sl]); m = len(ns)
    G = np.zeros(eng.nb); qlo = np.zeros(eng.nb); qhi = np.zeros(eng.nb)
    for j, b in enumerate(eng.reg_bus):
        G[b] = max(eng.qmax[j], -eng.qmin[j], 1.0) / slope; qlo[b] = eng.qmin[j]; qhi[b] = eng.qmax[j]
    def qreg(vm):
        raw = G * (eng.vset - vm); return np.clip(raw, qlo, qhi), (raw > qlo) & (raw < qhi)
    def F(V):
        S = V * np.conj(Y @ V); q, act = qreg(np.abs(V))
        return np.r_[(S.real - P)[ns], (S.imag - Qf - q)[ns]], act, q
    V = V0.copy(); f, act, q = F(V); err = np.abs(f).max()
    for it in range(maxit):
        if err < tol:
            return dict(ok=True, V=V, q=q, sat_hi=int(((q >= qhi - 1e-6) & (G > 0)).sum()), it=it)
        dth, dvm = eng._dSdV(Y, V)
        J = np.block([[dth[np.ix_(ns, ns)].real, dvm[np.ix_(ns, ns)].real],
                      [dth[np.ix_(ns, ns)].imag, dvm[np.ix_(ns, ns)].imag + np.diag((G * act)[ns])]])
        try: dx = np.linalg.solve(J, -f)
        except np.linalg.LinAlgError: return dict(ok=False)
        th0, vm0 = np.angle(V), np.abs(V); done = False
        # pas borné (0,03 pu et 0,1 rad au plus) : près du nez, le jacobien est presque singulier et un pas
        # entier fait sauter Newton vers une autre branche (tensions de 1,2 à 1,5 pu aux bus à condensateurs)
        a0 = min(1.0, 0.03 / max(np.abs(dx[m:]).max(), 1e-12), 0.1 / max(np.abs(dx[:m]).max(), 1e-12))
        for a in (a0, a0 / 2, a0 / 4, a0 / 8, a0 / 16, a0 / 32):
            th = th0.copy(); vm = vm0.copy(); th[ns] += a * dx[:m]; vm[ns] += a * dx[m:]
            Vt = vm * np.exp(1j * th); ft, actt, qt = F(Vt); e = np.abs(ft).max()
            if np.isfinite(e) and e < err: done = True; break
        if not done or vm.min() < 0.5 or vm.max() > 1.15: return dict(ok=False)
        V, f, act, q, err = Vt, ft, actt, qt, e
    return dict(ok=False)


def marge_statisme(eng, sn, k, c, slope=0.02, step=0.0025, lam_max=1.6, P0=None, Q0=None, lam_min=1.0):
    """Courbe P-V par continuation (départ à chaud) avec régulateurs à statisme. Retourne le DataFrame des pas réussis."""
    P0 = eng.P.loc[sn].values.copy() if P0 is None else P0; Q0 = eng.Q.loc[sn].values.copy() if Q0 is None else Q0
    Y = Ymat(eng, k, c)
    # départ : solution PV/PQ du moteur (même bassin que la solution d'exploitation), sinon départ DC
    saveP, saveQ = eng.P.loc[sn].copy(), eng.Q.loc[sn].copy()
    eng.P.loc[sn] = P0 * lam_min; eng.Q.loc[sn] = Q0 * lam_min
    r0 = eng.solve(sn, k=k, c=c); eng.P.loc[sn] = saveP; eng.Q.loc[sn] = saveQ
    if r0["ok"]:
        V = r0["V"]
    else:
        th = eng.dc_start(Y, P0); vm = np.ones(eng.nb); vm[eng.reg_bus] = eng.vset[eng.reg_bus]; vm[eng.sl] = eng.vset[eng.sl]
        V = vm * np.exp(1j * th)
    rows = []
    for lam in np.round(np.arange(lam_min, lam_max + 1e-9, step), 5):
        r = newton_statisme(eng, Y, P0 * lam, Q0 * lam, V, slope)
        if not r["ok"]: break
        V = r["V"]; v = np.abs(V); S = V * np.conj(Y @ V)
        rows.append(dict(lam=lam, v_min=v.min(), bus=eng.buses[v.argmin()], v_max=v.max(), sat_hi=r["sat_hi"], p_slack=S.real[eng.sl] - P0[eng.sl] * lam))
    return pd.DataFrame(rows)


def resume_marge(d, load, vlim=0.90):
    """lambda max (nez de la courbe) et lambda max avec toutes les tensions >= vlim."""
    if not len(d): return dict(lam_nez=np.nan, lam_090=np.nan, vmin_1=np.nan)
    ok = d[d.v_min >= vlim]; d1 = d[np.isclose(d.lam, 1.0)]
    return dict(lam_nez=float(d.lam.max()), lam_090=float(ok.lam.max()) if len(ok) else np.nan, vmin_1=float(d1.v_min.iloc[0]) if len(d1) else np.nan,
                GW_nez=round(load * d.lam.max() / 1e3, 1), sat_1=int(d1.sat_hi.iloc[0]) if len(d1) else -1, p_slack_1=float(d1.p_slack.iloc[0]) if len(d1) else np.nan)
