"""A2b : d'où vient l'écart entre l'équivalent de Kron calculé et celui présent dans le réseau final ?"""
from _commun import *
full = R.load(os.path.join(R.CURRENT, "elec_full.nc")); T = R.mercator_table(full)
src = R.load(os.path.join(R.CURRENT, "elec_solved.nc")); n, eng = final()
for reg, box in R.KRON_REGIONS.items():
    br, sh, info = R.kron_equivalent(full, src, T, box=box)
    rows = []
    for _, b in br.iterrows():
        nm = f"KRON {reg} {b.bus0}-{b.bus1}"
        if nm in n.lines.index:
            r = n.lines.loc[nm]; z = (r.r + 1j * r.x) / r.v_nom ** 2; kind = "Line"; extra = r.num_parallel
        else:
            r = n.transformers.loc[nm]; z = (r.r + 1j * r.x) / r.s_nom; kind = "Transfo"; extra = r.tap_ratio
        rows.append(dict(nom=nm, kind=kind, r_calc=b.z_pu.real, x_calc=b.z_pu.imag, r_final=z.real, x_final=z.imag, extra=extra))
    D = pd.DataFrame(rows); D["dx"] = D.x_final / D.x_calc - 1; D["r_ratio"] = D.r_final / D.r_calc.clip(lower=1e-6)
    print(reg, "| x changé (> 0,1 %) :", int((D.dx.abs() > 1e-3).sum()), "| r changé (> 1 %) :", int(((D.r_ratio - 1).abs() > 0.01).sum()), "sur", len(D))
    print(D[((D.r_ratio - 1).abs() > 0.01) | (D.dx.abs() > 1e-3)].groupby("kind").agg(n=("nom", "size"), r_ratio_med=("r_ratio", "median"), dx_med=("dx", "median")).to_string())
    print(D[D.kind == "Transfo"].head(4).to_string())
import inspect, run_pf
src_code = inspect.getsource(run_pf.prepare_for_ac_pf)
print([l.strip() for l in src_code.splitlines() if "transformers" in l][:12])
