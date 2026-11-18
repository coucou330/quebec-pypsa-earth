"""C3b : où le réseau complet casse-t-il ? Écarts d'angle par branche dans la solution DC (départ de Newton)."""
from complet import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full); src = R.load(os.path.join(R.CURRENT, "elec_solved.nc"))
_l0 = R.lopf
def lopf_ror(n, solver="highs", log=None):
    n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
    return _l0(n, solver=solver, log=log)
n, info = build_complet(full, T, lopf=lopf_ror)
eng = SparseEngine(n, info["ac"]["regulators"])
import pickle; pickle.dump((n, info), open("_complet_ror.pkl", "wb"))
Y = eng.Y0 + np.diag(1j * (0.37 * eng.b_react + eng.b_cap))
for nom, P in (("charges redescendues", eng.P.loc[PEAK].values),):
    th = np.degrees(eng.dc_start(Y, P)); ths = pd.Series(th, index=eng.buses)
    rows = []
    for kind, df in (("ligne", n.lines), ("transfo", n.transformers)):
        for nm, r in df.iterrows():
            if r.bus0 in eng.ix and r.bus1 in eng.ix:
                xpu = (r.x / r.v_nom ** 2) if kind == "ligne" else r.x / r.s_nom
                d = ths[r.bus0] - ths[r.bus1]
                rows.append(dict(nom=nm, kind=kind, bus0=r.bus0, bus1=r.bus1, kV=n.buses.v_nom[r.bus0], km=getattr(r, "length", 0.0), x_pu=xpu, dtheta=d, P_dc=np.radians(d) / xpu, Pmax=1 / xpu))
    D = pd.DataFrame(rows); D["abs"] = D.dtheta.abs()
    print(nom, ": étendue des angles DC", round(th.max() - th.min()), "° ; branches à plus de 30° :", int((D["abs"] > 30).sum()), "; plus de 60° :", int((D["abs"] > 60).sum()), "; plus de 90° :", int((D["abs"] > 90).sum()))
    print(D.sort_values("abs", ascending=False).head(25).round(3).drop(columns="abs").to_string(index=False))
    D.to_csv("c3b_angles_dc.csv", index=False)
    print("\npar tension, branches à plus de 30° :", D[D["abs"] > 30].groupby("kV").size().to_dict())
