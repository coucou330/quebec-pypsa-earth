"""A3 : cohérence état PV/PQ <-> V sur les 168 heures, Q et P du slack, total des inductances."""
from _commun import *
n, eng = final()
V = rd("tensions_pu.csv"); S = rd("regulateurs_etat.csv"); Qr = rd("regulateurs_Q_MVAr.csv"); Tt = rd("semaine_resume.csv")
K = rd("inductances_k.csv").reindex(columns=eng.buses).fillna(0); C = rd("condensateurs_c.csv").reindex(columns=eng.buses).fillna(0)
A = rd("angles_deg.csv")
vset = pd.Series(eng.vset, index=eng.buses); rb = {g: eng.buses[b] for g, b in zip(eng.regs, eng.reg_bus)}
qmax = pd.Series(eng.qmax, index=eng.regs); qmin = pd.Series(eng.qmin, index=eng.regs)
# 1. point fixe : résidu des équations avec V stocké et l'état stocké
res = []; comp = []; pvdev = []; qviol = []
for sn in V.index:
    Y = eng.Y0 + np.diag(1j * (K.loc[sn].values * eng.b_react + C.loc[sn].values * eng.b_cap))
    Vc = V.loc[sn, eng.buses].values * np.exp(1j * np.radians(A.loc[sn, eng.buses].values))
    Sc = Vc * np.conj(Y @ Vc)
    P = eng.P.loc[sn].values; Qf = eng.Q.loc[sn].values
    st = S.loc[sn, eng.regs].values.astype(int)
    mis_p = np.abs(np.delete(Sc.real - P, eng.sl)).max()
    Qs = Qf.copy(); pv = set()
    for j, b in enumerate(eng.reg_bus):
        if st[j] == 0: pv.add(b)
        else: Qs[b] += eng.qmax[j] if st[j] > 0 else eng.qmin[j]
    pq = [i for i in range(eng.nb) if i != eng.sl and i not in pv]
    mis_q = np.abs((Sc.imag - Qs)[pq]).max()
    res.append((sn, mis_p, mis_q))
    for j, b in enumerate(eng.reg_bus):
        v = abs(Vc[b]); q = Sc.imag[b] - Qf[b]
        if st[j] == 0:
            pvdev.append(abs(v - eng.vset[b]))
            if q > eng.qmax[j] + 1 or q < eng.qmin[j] - 1: qviol.append((sn, eng.regs[j], q))
        elif (st[j] > 0 and v > eng.vset[b] + 1e-3) or (st[j] < 0 and v < eng.vset[b] - 1e-3):
            comp.append((sn, eng.regs[j], st[j], round(v - eng.vset[b], 4)))
Rz = pd.DataFrame(res, columns=["sn", "mis_p", "mis_q"]).set_index("sn")
print(f"1. résidu max des équations avec (V, état) stockés : P {Rz.mis_p.max():.2e} MW, Q {Rz.mis_q.max():.2e} MVAr (tolérance du moteur 1e-3) ; écart max |V - consigne| des PV : {max(pvdev):.1e}")
print(f"2. régulateurs PV hors limites Q : {len(qviol)} ; violations de complémentarité (> 0,001 pu) : {len(comp)} -> {comp}")
# plusieurs régulateurs par bus ?
print("   régulateurs partageant un bus :", pd.Series(eng.reg_bus).duplicated().sum())
# 3. slack
g = n.generators; sl = g.index[g.control == "Slack"]; b = eng.slack_bus
print(f"3. slack : {list(sl)} au bus {b}, limites Q du slack {g.loc[sl, ['q_min', 'q_max']].values.tolist()}")
print(f"   Q slack sur la semaine : {Tt.q_slack.min():.0f} à {Tt.q_slack.max():.0f} MVAr ; heures hors limites : {int(((Tt.q_slack > g.loc[sl, 'q_max'].sum()) | (Tt.q_slack < g.loc[sl, 'q_min'].sum())).sum())}")
gens_b = g[(g.bus == b) & ~g.carrier.isin(["regulator", "load_shedding"])]; su_b = n.storage_units[n.storage_units.bus == b]
pcap = su_b.p_nom.sum() + gens_b[~gens_b.index.str.contains("placeholder")].p_nom.sum()
# P du LOPF aux unités du bus slack : déjà dans eng.P (sauf le générateur slack lui-même) ; p_slack = ce que le slack ajoute
su_p = n.storage_units_t.p_set.reindex(columns=su_b.index).fillna(0).sum(axis=1).reindex(Tt.index)
tot = Tt.p_slack + su_p
load = n.loads_t.p_set.sum(axis=1)
print(f"   unités réelles au bus slack : {list(su_b.index)}, {pcap:.0f} MW installés ; dispatch LOPF de ces unités : {su_p.min():.0f} à {su_p.max():.0f} MW")
print(f"   P ajoutée par le slack (= pertes AC, le LOPF est sans pertes) : {Tt.p_slack.min():.0f} à {Tt.p_slack.max():.0f} MW ; à la pointe {Tt.p_slack[PEAK]:.0f} MW = {Tt.p_slack[PEAK] / load[PEAK]:.1%} de la charge")
print(f"   P totale produite au bus slack : max {tot.max():.0f} MW pour {pcap:.0f} MW installés ; heures au-dessus de l'installé : {int((tot > pcap + 1).sum())}/168 ; à la pointe {tot[PEAK]:.0f} MW")
# 4. inductances : total à k = 1 et enclenché, contre HQ (25 500 Mvar installés, ~11 000 à la pointe)
tot_k1 = -(eng.b_react).sum() / 1e3
print(f"4. inductances : {tot_k1:.1f} GVAr installés à k = 1 dans le modèle (HQ 1996 : 25,5 GVAr) ; enclenché : min {Tt.q_react_MVAr.min() / 1e3:.1f}, max {Tt.q_react_MVAr.max() / 1e3:.1f} GVAr ; "
      f"heures où l'enclenché dépasse 25,5 GVAr : {int((Tt.q_react_MVAr > 25500).sum())}/168")
L = n.lines; bq = (L.b * L.v_nom ** 2); print(f"   charge capacitive des lignes à 1 pu : 735 kV {bq[L.v_nom >= 700].sum() / 1e3:.1f} GVAr, 315 kV {bq[(L.v_nom < 700)].sum() / 1e3:.1f} GVAr ; par km de 735 kV : {bq[L.v_nom >= 700].sum() / (L.length * L.num_parallel)[L.v_nom >= 700].sum():.2f} MVAr/km")
print(f"   condensateurs : {eng.b_cap.sum() / 1e3:.1f} GVAr installés ; enclenchés à la pointe {Tt.q_cap_MVAr[PEAK] / 1e3:.1f} ; au creux {Tt.q_cap_MVAr.min() / 1e3:.1f}")
