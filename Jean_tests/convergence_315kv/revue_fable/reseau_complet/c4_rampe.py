"""C4 : rampe de charge sur le réseau complet (toutes les injections × lambda, condensateurs c = lambda,
réactif des régulateurs illimité) : jusqu'où Newton converge-t-il, et quels bus lâchent ?"""
from complet import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
ratio = float(sys.argv[1]) if len(sys.argv) > 1 else 0.10
n, info = build_complet(full, T, lopf=lopf_ror, max_ratio=ratio)
eng = SparseEngine(n, info["ac"]["regulators"])
P1, Q1 = eng.P.loc[PEAK].values.copy(), eng.Q.loc[PEAK].values.copy()
eng.qmin[:] = -1e12; eng.qmax[:] = 1e12
Vs = None; vm_last = None
for lam in [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.85, 0.9, 0.95, 1.0]:
    eng.P.loc[PEAK] = P1 * lam; eng.Q.loc[PEAK] = Q1 * lam
    r = eng.solve(PEAK, k=0.6 - 0.3 * lam, c=lam, Vstart=Vs)
    if not r["ok"]: r = eng.solve(PEAK, k=0.6 - 0.3 * lam, c=lam)
    if not r["ok"]: print(f"lambda = {lam} : Newton diverge"); break
    vm = pd.Series(np.abs(r["V"]), index=eng.buses); Vs = r["V"]; vm_last = vm
    lo = vm.sort_values().head(4)
    print(f"lambda = {lam} : V {vm.min():.3f}-{vm.max():.3f} ; bus < 0,90 : {(vm < 0.9).sum()} ; > 1,10 : {(vm > 1.1).sum()} ; plus bas : {[(b, int(n.buses.v_nom[b]), round(v, 2)) for b, v in lo.items()]} ; plus hauts : {[(b, int(n.buses.v_nom[b]), round(v, 2)) for b, v in vm.sort_values().tail(3).items()]}", flush=True)
eng.P.loc[PEAK] = P1; eng.Q.loc[PEAK] = Q1
if vm_last is not None:
    vm_last.to_csv("c4_tensions_dernier_pas.csv")
    Sload = n.loads_t.p_set.loc[PEAK].groupby(n.loads.bus).sum()
    lo = vm_last.sort_values().head(15)
    print(pd.DataFrame(dict(kV=n.buses.v_nom[lo.index], V=lo.round(3), charge_MW=Sload.reindex(lo.index).fillna(0).round(), lat=n.buses.y[lo.index].round(2), lon=n.buses.x[lo.index].round(2))).to_string())
