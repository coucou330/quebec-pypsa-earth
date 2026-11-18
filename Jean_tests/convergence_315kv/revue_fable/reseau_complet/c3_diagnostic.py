"""C3 : pourquoi le réseau complet diverge-t-il à la pointe ?
On fait redescendre les charges progressivement : phi = 0 (toutes les charges aux bus >= 315 kV, comme le
modèle 315 kV, le réseau < 315 kV n'étant qu'un réseau passif en parallèle) -> phi = 1 (toutes les charges
à leur bus d'origine). On cherche le phi maximal où Newton converge et les bus qui lâchent."""
from complet import *
import time, pickle
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full); src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
n, info = build_complet(full, T, lopf=lopf_ror)
eng = SparseEngine(n, info["ac"]["regulators"]); nb = eng.nb
up_bus = src.loads.bus.reindex(n.loads.index); dn_bus = n.loads.bus
lp = n.loads_t.p_set.loc[PEAK].reindex(n.loads.index); lq = n.loads_t.q_set.loc[PEAK].reindex(n.loads.index)
dP = np.zeros(nb); dQ = np.zeros(nb)     # ce qu'il faut ajouter à l'injection « descendue » pour revenir à l'état « remonté »
for l in n.loads.index:
    a, b = dn_bus[l], up_bus[l]
    if a != b and a in eng.ix and b in eng.ix:
        dP[eng.ix[a]] += lp[l]; dP[eng.ix[b]] -= lp[l]; dQ[eng.ix[a]] += lq[l]; dQ[eng.ix[b]] -= lq[l]
P1, Q1 = eng.P.loc[PEAK].values.copy(), eng.Q.loc[PEAK].values.copy()
qmin0, qmax0 = eng.qmin.copy(), eng.qmax.copy(); eng.qmin[:] = -1e12; eng.qmax[:] = 1e12   # réactif illimité : on teste le réseau seul
def essai(phi, k=0.37, c=1.0, Vs=None):
    eng.P.loc[PEAK] = P1 + (1 - phi) * dP; eng.Q.loc[PEAK] = Q1 + (1 - phi) * dQ
    r = eng.solve(PEAK, k=k, c=c, Vstart=Vs); eng.P.loc[PEAK] = P1; eng.Q.loc[PEAK] = Q1; return r
Vs = None; last = None
for phi in [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0]:
    r = essai(phi, Vs=Vs)
    if not r["ok"]: r = essai(phi)
    if not r["ok"]: print(f"phi = {phi:.1f} : Newton diverge"); break
    vm = np.abs(r["V"]); Vs = r["V"]; last = (phi, vm.copy())
    print(f"phi = {phi:.1f} : converge ; V {vm.min():.3f} (bus {eng.buses[vm.argmin()]}, {eng.vnom[vm.argmin()]:.0f} kV) à {vm.max():.3f} ; bus < 0,90 : {(vm < 0.9).sum()} ; < 0,80 : {(vm < 0.8).sum()} ; pertes {r['p_slack']:.0f} MW", flush=True)
eng.qmin[:], eng.qmax[:] = qmin0, qmax0
# bus de charge : charge à la pointe contre « force » du réseau local (somme des V^2/x des branches du bus)
B = n.buses; Sload = n.loads_t.p_set.loc[PEAK].groupby(n.loads.bus).sum()
cap = pd.Series(0.0, index=B.index)
for df, zf in ((n.lines, lambda r: r.x / r.v_nom ** 2), (n.transformers, lambda r: r.x / r.s_nom)):
    for _, r in df.iterrows():
        y = 1 / abs(zf(r)) if zf(r) else 0
        cap[r.bus0] += y; cap[r.bus1] += y
D = pd.DataFrame(dict(v_nom=B.v_nom, charge_MW=Sload.reindex(B.index).fillna(0), force_MW=cap, degre=pd.concat([n.lines.bus0, n.lines.bus1, n.transformers.bus0, n.transformers.bus1]).value_counts().reindex(B.index).fillna(0)))
D = D[(D.v_nom < 315) & D.index.isin(eng.buses)]; D["charge_sur_force"] = D.charge_MW / D.force_MW
if last is not None: D["V_dernier_phi"] = pd.Series(last[1], index=eng.buses).reindex(D.index)
D.sort_values("charge_sur_force", ascending=False).to_csv("c3_bus_charge.csv")
print("\ncharge à la pointe par tension de bus (MW) :", D.groupby("v_nom").charge_MW.sum().round().to_dict())
print("charge par bus < 315 kV : médiane %.0f MW, max %.0f MW ; bus de plus de 300 MW : %d ; de plus de 500 MW : %d" % (D.charge_MW[D.charge_MW > 0].median(), D.charge_MW.max(), (D.charge_MW > 300).sum(), (D.charge_MW > 500).sum()))
print("bus dont la charge dépasse 30 % de V^2/x de leurs branches :", int((D.charge_sur_force > 0.3).sum()), "; 50 % :", int((D.charge_sur_force > 0.5).sum()), "; antennes (1 seule branche) avec charge :", int(((D.degre <= 1) & (D.charge_MW > 0)).sum()))
print(D.sort_values("charge_sur_force", ascending=False).head(15).round(3).to_string())
# transformateurs 315/735 -> < 315 kV : combien, et quelle charge passe dessous ?
tr = n.transformers[n.transformers.index.str.startswith("BT ")]
hv = tr[(B.v_nom[tr.bus0].values >= 315) | (B.v_nom[tr.bus1].values >= 315)]
print(f"\ntransformateurs entre le réseau >= 315 kV et le réseau < 315 kV : {len(hv)} (2 000 MVA chacun) pour {D.charge_MW.sum():.0f} MW de charge dessous ; postes >= 315 kV concernés : {len(set(hv.bus0) | set(hv.bus1)) - len(set(hv.bus0) & set(B.index[B.v_nom < 315])) - len(set(hv.bus1) & set(B.index[B.v_nom < 315]))}")
