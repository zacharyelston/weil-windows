"""Compact data for the 'window minimum' section of the zero-torus page (docs/viz/zero_torus.html).

Reads committed R2 results only:
  data/research_r2_rigorous/kb_cert.json, hermite_cert.json, summary.json
Writes data/viz/r2_panel.json. rh_zero_torus substitutes it for __R2__ in the page template.

Per object and sector (degree 1, c = 2*pi*x/q):
  kb    certified Kaiser-Bessel upper bound (Theorem C; better of beta = c, c - 2)   [C]
  herm  certified Hermite upper bound (Theorem C)                                  [C]
  ritz  finite-basis minimum from rh2's scans (an upper bound, numerical)          [N]
  norm  ln ||f_win||^2 - 2*beta for the beta = c - 2 trial (certified lower bound;  [C]
        GPT's upper enclosure agrees to 1.5e-4 relative), with the power-law fit   [N]
ζ supports: certified lower bound (certificates), KB and Hermite bounds, KB trial RQ, and the KB cutoff
  Tc = A/B (beta = c - 2): Theorem 2 charges each zero min(B, A/|gamma|)^2, i.e. B^2 below Tc.
"""
import json, math, os
from flint import arb

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
D = os.path.join(ROOT, "data", "research_r2_rigorous")


def mid(s):
    s = str(s)
    return float(arb(s.split("+/-")[0].strip(" ["))) if "+/-" in s else float(s)


def load(name):
    with open(os.path.join(D, name)) as fh:
        return json.load(fh)["rows"]


def fnum(v):
    return None if v in (None, "None") else float(v)


def lsq2(xs, ys):
    n = len(xs); sx = sum(xs); sy = sum(ys); sxx = sum(x * x for x in xs); sxy = sum(x * y for x, y in zip(xs, ys))
    b = (n * sxy - sx * sy) / (n * sxx - sx * sx); a = (sy - b * sx) / n
    return a, b


grid = {}
for r in load("kb_cert.json"):
    D_, s, xq = int(r["D"]), int(r["s"]), int(r["xq"])
    q = int(r["q"]); c = 2 * math.pi * mid(r["x"]) / q
    cand = [k for k in r["candidates"] if int(k["dbeta"]) == 2 and int(k["m"]) == 2][0]
    beta = mid(cand["beta"])
    grid[(D_, s, xq)] = {"c": round(c, 4), "kb": float(r["bound_upper"]), "ritz": fnum(r["measured_ritz"]),
                         "norm": round(math.log(mid(cand["norm2"])) - 2 * beta, 5),
                         "a": round(mid(cand["a_t2"]), 6) if D_ == 1 else None}
herm = {}
for r in load("hermite_cert.json"):
    D_, s, xq = int(r["D"]), int(r["s"]), int(r["xq"])
    q = int(r["q"]); c = 2 * math.pi * mid(r["x"]) / q
    herm[(D_, s, xq)] = {"c": round(c, 4), "herm": float(r["bound_upper"]), "ritz": fnum(r["measured_ritz"])}

objects = [1, 5, 8, -3, -4, -7, -20]
out = {"objects": [], "supports": [], "source": "data/research_r2_rigorous (Theorem C), docs/AUDIT_R2_CONJECTURE_R.md"}
for D_ in objects:
    for s in (0, 1):
        rows = []
        xqs = sorted({k[2] for k in herm if k[0] == D_ and k[1] == s} | {k[2] for k in grid if k[0] == D_ and k[1] == s})
        for xq in xqs:
            h = herm.get((D_, s, xq)); g = grid.get((D_, s, xq))
            ritz = (g or {}).get("ritz") or (h or {}).get("ritz")
            row = {"xq": xq, "c": (g or h)["c"]}
            if ritz is not None: row["ritz"] = ritz
            if h: row["herm"] = h["herm"]
            if g: row.update(kb=g["kb"], norm=g["norm"])
            if g and g["a"] is not None: row["a"] = g["a"]
            rows.append(row)
        nr = [(math.log(r["c"]), r["norm"]) for r in rows if "norm" in r]
        a0, b0 = lsq2([x for x, _ in nr], [y for _, y in nr])
        resid = max(abs(y - a0 - b0 * x) for x, y in nr)
        out["objects"].append({"D": D_, "s": s, "parity": "even" if s == 0 else "odd", "rows": rows,
                               "fit": {"a": round(a0, 4), "p2": round(-b0, 3), "maxres": round(resid, 3)}})
tc = {}
for r in load("kb_cert_supports.json"):
    cand = [k for k in r["candidates"] if int(k["dbeta"]) == 2 and int(k["m"]) == 2][0]
    tc[(str(r["support_2a"]), r["parity"])] = round(mid(cand["Tc"]), 3)  # A/B: full charge B^2 below, (Tc/gamma)^2 B^2 above
for r in json.load(open(os.path.join(D, "summary.json")))["supports"]:
    out["supports"].append({"support": r["support_2a"], "x": round(r["x"], 4), "Tc": tc[(str(r["support_2a"]), r["parity"])],
                            "c": round(2 * math.pi * r["x"], 4), "parity": r["parity"],
                            "cert_lower": float(r["cert_lower"]), "kb": float(r["kb_B_cert"]),
                            "herm": float(r["hermite_B_cert"]), "kb_rq": float(r["kb_rq"])})
os.makedirs(os.path.join(ROOT, "data", "viz"), exist_ok=True)
with open(os.path.join(ROOT, "data", "viz", "r2_panel.json"), "w") as fh:
    json.dump(out, fh, separators=(",", ":"))
print("wrote data/viz/r2_panel.json:", len(out["objects"]), "classes,", len(out["supports"]), "support rows")
for o in out["objects"]:
    print(o["D"], o["parity"], o["fit"], len(o["rows"]))
