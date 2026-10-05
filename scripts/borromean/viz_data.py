"""Data for the 'Borromean primes' section of the zero-torus page, from data/borromean/step1.json (committed results).
Writes data/viz/borromean_panel.json; rh-zero-torus substitutes it for __BORROMEAN__ in the page template."""
import json, os
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
d = json.load(open(os.path.join(ROOT, "data", "borromean", "step1.json")))
g3 = d["g3"]
neg = {r["p"]: r["R_B"] for r in g3["negative_control"]}
rows = [{"p": r["p"], "a": r["a_p"], "R": round(r["R"], 4), "borromean": r["borromean"],
         **({"RB": round(neg[r["p"]], 4)} if r["p"] in neg else {})} for r in g3["rows"]]
out = {"pair": d["pair"], "T": d["T"], "n_zeros_rho": d["rho"]["g1_38"]["n_zeros"], "rows": rows,
       "max_dev_unlinked": round(g3["max_dev_unlinked"], 4), "field": d["g0"]["field"], "conductor": d["g0"]["conductor"],
       "source": "data/borromean/step1.json (issue the private research log (issue 7), step 1)"}
os.makedirs(os.path.join(ROOT, "data", "viz"), exist_ok=True)
json.dump(out, open(os.path.join(ROOT, "data", "viz", "borromean_panel.json"), "w"), separators=(",", ":"))
print("rows", len(rows), "borromean", [r["p"] for r in rows if r["borromean"]], "a=-2 non-Borromean", [r["p"] for r in rows if r["a"] == -2 and not r["borromean"]][:8])
