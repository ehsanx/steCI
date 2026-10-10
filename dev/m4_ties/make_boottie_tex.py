"""Fill dev/m4_ties/boottie_appendix_template.tex from se_subset outputs -> boottie_appendix.tex.

    python dev/m4_ties/make_boottie_tex.py

Table rows come from se_subset_table.tex (written by se_subset.R); the Table-3 SE-ratio
figures come from se_subset.log's 'implied Table e1' block. The output is the exact block
pasted into appendix/appendix.tex as the last subsection of Web Appendix D.
"""
import pathlib, re

HERE = pathlib.Path(__file__).resolve().parent
tpl = (HERE / "boottie_appendix_template.tex").read_text(encoding="utf-8")
rows = (HERE / "se_subset_table.tex").read_text(encoding="utf-8").splitlines()
ia, ib = rows.index("% Panel A rows (production, balanced)"), rows.index("% Panel B rows (pilot)")
panel_a, panel_b = rows[ia + 1:ib], rows[ib + 1:]
assert len(panel_a) == 18 and len(panel_b) == 6, (len(panel_a), len(panel_b))

log = (HERE / "se_subset.log").read_text(encoding="utf-8-sig")
blk = log[log.index("implied Table e1"):]
vals = {}
for n in ("300", "1000"):
    m = re.search(r"^\s*%s\s+27\s+\d+\s+([0-9.]+)\s+(-?[0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)\s+([0-9.]+)" % n, blk, re.M)
    assert m, n
    pub, shift, mcse, implied, pnow, pnew = (float(g) for g in m.groups())
    assert abs(pub + shift - implied) < 1e-3
    vals.update({f"SER{n}": f"{implied:.3f}", f"PUB{n}": f"{pub:.3f}", f"MC{n}": f"{mcse:.3f}",
                 f"P{n}": f"{pnow:.2f}", f"Q{n}": f"{pnew:.2f}"})

out = tpl.replace("@PANELA@", "\n".join(panel_a)).replace("@PANELB@", "\n".join(panel_b))
for k, v in vals.items():
    out = out.replace(f"@{k}@", v)
assert not re.findall(r"@[A-Z0-9]+@", out), re.findall(r"@[A-Z0-9]+@", out)
(HERE / "boottie_appendix.tex").write_text(out, encoding="utf-8")
print("wrote boottie_appendix.tex;", vals)
