"""C2 : le réseau complet converge-t-il à la pointe ? Avec quel profil de tension ? Quelle marge ?
Deux dispatchs : celui du cas final (LOPF tel quel) et celui avec fil de l'eau prioritaire (constat B4)."""
from complet import *
import time
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
which = sys.argv[1] if len(sys.argv) > 1 else "ror"
ratio = float(sys.argv[2]) if len(sys.argv) > 2 else 0.10
t = time.time(); n, info = build_complet(full, T, lopf=lopf_ror if which == "ror" else None, max_ratio=ratio)
print(f"[{which}] construit en {time.time() - t:.0f} s :", info["complet"], "| délestage LOPF %.4f %%" % info["lopf"]["shed_pct"], flush=True)
eng = SparseEngine(n, info["ac"]["regulators"])
print("bus du réseau AC principal :", eng.nb, "| par tension :", pd.Series(eng.vnom).value_counts().sort_index().to_dict())
print("sous-réseaux :", n.buses.sub_network.value_counts().head(6).to_dict())
load = n.loads_t.p_set.sum(axis=1)
# 1. Newton seul, régulateurs sans limites, à plusieurs états de shunts
qmin0, qmax0 = eng.qmin.copy(), eng.qmax.copy()
for k, c in [(0.37, 1.0), (0.2, 1.0), (0.37, 0.5), (0.5, 0.0)]:
    eng.qmin[:] = -1e12; eng.qmax[:] = 1e12
    t = time.time(); r = eng.solve(PEAK, k=k, c=c); eng.qmin[:], eng.qmax[:] = qmin0, qmax0
    if r["ok"]:
        vm = np.abs(r["V"]); over = int(((r["q_reg"] > qmax0 + 1) | (r["q_reg"] < qmin0 - 1)).sum())
        print(f"pointe, Q illimité, k={k}, c={c} : Newton converge ; V {vm.min():.3f}-{vm.max():.3f} ; bus < 0,90 : {(vm < 0.9).sum()}, > 1,10 : {(vm > 1.1).sum()} ; régulateurs au-delà de leurs limites : {over}/{len(eng.regs)} ; pertes {r['p_slack']:.0f} MW ({time.time() - t:.1f} s)", flush=True)
    else: print(f"pointe, Q illimité, k={k}, c={c} : Newton diverge", flush=True)
# 2. avec limites Q et manœuvre des shunts (même logique que le 315 kV)
t = time.time(); r = R.solve_hour_reserve(eng, PEAK, R.K_PEAK, 1.0)
if r["ok"]:
    tab, vm = bilan(n, eng, r)
    print(f"\npointe, limites Q + shunts manœuvrés : Newton converge ; régulateurs PV hors limites : {r['n_viol']} ; en butée : {int((r['state'] != 0).sum())}/{len(eng.regs)} ; V {vm.min():.3f}-{vm.max():.3f} ; "
          f"bus hors 0,90-1,10 : {int(((vm < 0.9) | (vm > 1.1)).sum())}/{eng.nb} ; hors 0,95-1,05 : {int(((vm < 0.95) | (vm > 1.05)).sum())} ; pertes {r['p_slack']:.0f} MW ({time.time() - t:.0f} s)")
    print(tab.to_string()); 
    low = vm.sort_values().head(12); print("bus les plus bas :", [(b, n.buses.v_nom[b], round(v, 3)) for b, v in low.items()])
    vm.to_csv(f"c2_tensions_pointe_{which}.csv")
    import pickle; pickle.dump(dict(k=r["k"], c=r["c"], V=r["V"], state=r["state"]), open(f"_c2_{which}.pkl", "wb"))
    # 3. marge de charge (critère : Newton + limites Q ; shunts re-manœuvrés)
    t = time.time(); d = R.load_margin(eng, PEAK, mode="remanoeuvre", step=0.01, compl_tol=9.0)
    print(f"marge (Newton + limites Q) : lambda max {d.lam.max() if len(d) else float('nan')} ; avec V min >= 0,90 : {d.lam[d.v_min >= 0.9].max() if len(d) and (d.v_min >= 0.9).any() else float('nan')} ({time.time() - t:.0f} s)")
    print(d.round(3).to_string(index=False))
else:
    print("\npointe, limites Q + shunts manœuvrés : pas de solution", r.get("iters"))
