"""C5 : jusqu'où peut-on descendre ? Élagage automatique.
On monte la charge par paliers de 5 % (réactif des régulateurs illimité). Quand Newton diverge, les charges des
bus < 315 kV dont la tension s'effondre dans le palier qui diverge (V < 0,90 après 6 itérations de Newton
depuis le dernier palier convergé ; à défaut les 5 plus bas) sont laissées à leur bus
>= 315 kV (comme dans la réduction de Nathan), et on recommence. On s'arrête quand la pointe converge.
Écrit c5_bus_exclus.csv (bus, tour, charge, tension) et c5_journal.txt."""
from complet import *
import time, pickle
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
which = sys.argv[1] if len(sys.argv) > 1 else "ror"
excl = {}; journal = open(f"c5_journal_{which}.txt", "w")
def say(*a):
    s = " ".join(str(x) for x in a); print(s, flush=True); journal.write(s + "\n"); journal.flush()
RAMPE = list(np.round(np.arange(0.2, 1.0001, 0.05), 2))
for tour in range(40):
    n, info = build_complet(full, T, lopf=lopf_ror if which == "ror" else None, exclure=set(excl))
    eng = SparseEngine(n, info["ac"]["regulators"])
    P1, Q1 = eng.P.loc[PEAK].values.copy(), eng.Q.loc[PEAK].values.copy(); q0 = (eng.qmin.copy(), eng.qmax.copy())
    eng.qmin[:] = -1e12; eng.qmax[:] = 1e12
    Vs = None; vm = None; ok_all = True; lam_ok = 0.0
    for lam in RAMPE:
        eng.P.loc[PEAK] = P1 * lam; eng.Q.loc[PEAK] = Q1 * lam
        r = eng.solve(PEAK, k=0.6 - 0.25 * lam, c=lam, Vstart=Vs)
        if not r["ok"]: r = eng.solve(PEAK, k=0.6 - 0.25 * lam, c=lam)
        if not r["ok"]:
            ok_all = False
            # autopsie du palier qui diverge : quelques itérations de Newton depuis le dernier palier convergé ;
            # les bus dont la tension s'effondre désignent la poche fautive
            if Vs is not None:
                kk = 0.6 - 0.25 * lam; Y = eng.Y0 + np.diag(1j * (kk * eng.b_react + lam * eng.b_cap))
                pv = sorted(set(eng.reg_bus)); pvs = set(pv); pq = [i for i in range(eng.nb) if i != eng.sl and i not in pvs]
                Vt, _, _ = eng.newton(Y, P1 * lam + 1j * Q1 * lam, Vs, np.array(pv), np.array(pq), maxit=6)
                vm = pd.Series(np.abs(Vt), index=eng.buses)
            break
        Vs = r["V"]; vm = pd.Series(np.abs(Vs), index=eng.buses); lam_ok = lam
    eng.P.loc[PEAK] = P1; eng.Q.loc[PEAK] = Q1; eng.qmin[:], eng.qmax[:] = q0
    c = info["complet"]
    say(f"tour {tour} : {len(excl)} bus exclus ; {c['MW_redescendus_pointe']:.0f} MW redescendus, {c['MW_restes']:.0f} MW restés au 315 kV ; Newton converge jusqu'à lambda = {lam_ok}"
        + (f" ; V {vm.min():.3f}-{vm.max():.3f}" if vm is not None else ""))
    if ok_all: break
    Sload = n.loads_t.p_set.loc[PEAK].groupby(n.loads.bus).sum()
    if vm is None:
        cand = Sload[[b for b in Sload.index if n.buses.v_nom[b] < 315]].sort_values(ascending=False).head(10).index
        low = pd.Series(np.nan, index=cand)
    else:
        lowv = vm[[b for b in vm.index if n.buses.v_nom[b] < 315 and Sload.get(b, 0) > 0]].sort_values()
        low = lowv[lowv < 0.90] if (lowv < 0.90).any() else lowv.head(5)
        low = low.head(25)
    for b, v in low.items(): excl[b] = dict(tour=tour, kV=n.buses.v_nom[b], charge_MW=round(float(Sload.get(b, 0))), V=v, lat=n.buses.y[b], lon=n.buses.x[b])
    pd.DataFrame(excl).T.to_csv(f"c5_bus_exclus_{which}.csv")
pickle.dump(sorted(excl), open(f"_c5_exclus_{which}.pkl", "wb"))
say("fin :", "la pointe converge (réactif illimité)" if ok_all else "pas de convergence après 40 tours")
journal.close()
