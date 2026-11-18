"""C4b : rampe fine entre 20 % et 40 % de la pointe, et autopsie du pas qui diverge (plus gros résidus)."""
from complet import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
n, info = build_complet(full, T, lopf=lopf_ror)
eng = SparseEngine(n, info["ac"]["regulators"])
P1, Q1 = eng.P.loc[PEAK].values.copy(), eng.Q.loc[PEAK].values.copy(); eng.qmin[:] = -1e12; eng.qmax[:] = 1e12
Vs = None; B = n.buses
for lam in np.round(np.arange(0.20, 0.42, 0.02), 2):
    k = 0.6 - 0.25 * lam; c = lam
    eng.P.loc[PEAK] = P1 * lam; eng.Q.loc[PEAK] = Q1 * lam
    r = eng.solve(PEAK, k=k, c=c, Vstart=Vs)
    if r["ok"]:
        Vs = r["V"]; vm = pd.Series(np.abs(Vs), index=eng.buses); ang = np.degrees(np.angle(Vs))
        print(f"lambda {lam} ok : V {vm.min():.3f} ({vm.idxmin()}) - {vm.max():.3f} ({vm.idxmax()}) ; angles {ang.min():.0f} à {ang.max():.0f}°", flush=True); continue
    # autopsie : un Newton depuis la dernière solution, on regarde les résidus après quelques itérations
    Y = eng.Y0 + np.diag(1j * (k * eng.b_react + c * eng.b_cap)); S = P1 * lam + 1j * Q1 * lam
    pv = sorted(set(eng.reg_bus)); pq = [i for i in range(eng.nb) if i != eng.sl and i not in set(pv)]
    V = Vs.copy(); Vm = np.abs(V); Vm[pv] = eng.vset[pv]; V = Vm * np.exp(1j * np.angle(V))
    for mi in (2, 4, 8, 30):
        Vt, ok, it = eng.newton(Y, S, V, np.array(pv), np.array(pq), maxit=mi)
        mis = Vt * np.conj(Y @ Vt) - S; mp = pd.Series(np.abs(mis.real), index=eng.buses); mq = pd.Series(np.abs(mis.imag), index=eng.buses)
        mq.loc[[eng.buses[i] for i in pv] + [eng.slack_bus]] = 0; mp.loc[eng.slack_bus] = 0
        vm = pd.Series(np.abs(Vt), index=eng.buses)
        print(f"lambda {lam} : après {mi} it., ok={ok} ; résidu P max {mp.max():.1f} au bus {mp.idxmax()} ({B.v_nom[mp.idxmax()]:.0f} kV), Q max {mq.max():.1f} au bus {mq.idxmax()} ({B.v_nom[mq.idxmax()]:.0f} kV) ; V {vm.min():.2f} ({vm.idxmin()}) - {vm.max():.2f} ({vm.idxmax()})")
    worst = (mp + mq).sort_values(ascending=False).head(8)
    Sload = n.loads_t.p_set.loc[PEAK].groupby(n.loads.bus).sum()
    print(pd.DataFrame(dict(kV=B.v_nom[worst.index], residu=worst.round(1), V=vm[worst.index].round(3), charge_pointe=Sload.reindex(worst.index).fillna(0).round(), lat=B.y[worst.index].round(2), lon=B.x[worst.index].round(2),
                            reg=[("REG " + b) in eng.regs for b in worst.index])).to_string())
    break
