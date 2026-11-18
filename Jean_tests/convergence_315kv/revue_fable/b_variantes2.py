"""Partie B, suite : variantes ajoutées après le constat B4 (le LOPF laisse 3,4 GW de fil de l'eau du sud à
l'arrêt à la pointe, parce que toute l'hydraulique a le même coût nul). « Fil de l'eau prioritaire » = coût
marginal de -0,2 $/MWh pour les centrales au fil de l'eau (elles ne peuvent pas stocker l'eau), comme le
modèle le fait déjà pour l'éolien (-0,1). Écrit b_variantes_part9.csv."""
import sys, time
import b_variantes as BV
from b_variantes import *

src = R.load(os.path.join(R.CURRENT, "elec_solved.nc")); BV.eng_ref_load = src.loads_t.p_set.sum(axis=1).reindex(BV.REF.index)
FAC = 1 - (BV.REF.p_slack / BV.eng_ref_load).clip(0, 0.1)

def lopf_combo(ror=True, charges=False):
    def f(n, solver="highs", log=None):
        if ror: n.generators.loc[n.generators.carrier == "ror", "marginal_cost"] = -0.2
        if charges: n.loads_t.p_set = n.loads_t.p_set.mul(FAC.reindex(n.snapshots).values, axis=0)
        return BV._lopf0(n, solver=solver, log=log)
    return f

QC = {"Québec": dict(lat=(44.0, 63.0), lon=(-80.0, -56.0))}
VAR2 = [
    ("B4 fil de l'eau prioritaire dans le LOPF", dict(lopf=lopf_combo()), None),
    ("B4 + B5 : fil de l'eau prioritaire, pertes réparties", dict(lopf=lopf_combo()), "pertes"),
    ("B4 + B5 : fil de l'eau prioritaire, charges réduites des pertes, pertes réparties", dict(lopf=lopf_combo(charges=True)), "pertes"),
    ("B4 + B5 + A1b : idem, sans la correction des circuits E3", dict(lopf=lopf_combo(charges=True), build=dict(circuits=False)), "pertes"),
    ("B4 + B5 + B2 : idem (avec E3), Kron tout le Québec + règle hybride", dict(lopf=lopf_combo(charges=True), regions=QC, T="hyb"), "pertes"),
]
if __name__ == "__main__":
    out = "b_variantes_part9.csv"; done = pd.read_csv(out).variante.tolist() if os.path.exists(out) else []
    for label, opt, post in VAR2:
        if label in done: continue
        try: row = BV.run(label, opt, post)
        except Exception as e:
            import traceback; row = dict(variante=label, erreur=repr(e)[:300]); traceback.print_exc()
        print(row, flush=True)
        pd.DataFrame([row]).to_csv(out, mode="a", header=not os.path.exists(out), index=False)
