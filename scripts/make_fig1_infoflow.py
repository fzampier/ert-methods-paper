#!/usr/bin/env python3
"""Figure 1 (information flow for one e-RT update), rebuilt as a script for round 2.
Same four-step panel as the author-drawn first version; wording follows the article:
step 2 is 'Wager recorded' (no 'locked'), no dashes, 'trial sites'. Output: figures/fig1_infoflow.{png,pdf}"""
import pathlib, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
plt.rcParams.update({"font.family": "Times New Roman", "mathtext.fontset": "stix", "pdf.fonttype": 42})
out = pathlib.Path(__file__).resolve().parents[1] / "figures"; out.mkdir(exist_ok=True)
INK, GREY, RULE, SOFT = "#111111", "#6b6b6b", "#9a9a9a", "#d9d9d9"
fig = plt.figure(figsize=(6.6, 6.7)); ax = fig.add_axes([0, 0, 1, 1]); ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.axis("off")
ax.add_patch(Rectangle((0.004, 0.004), 0.992, 0.992, fill=False, ec="#c8c8c8", lw=1.2))
XN, XT, XR, XE = 0.065, 0.165, 0.718, 0.935
def rules(y, color, lw):
    for x0, x1 in ((XN, 0.133), (XT, 0.686), (XR, XE)): ax.plot([x0, x1], [y, y], color=color, lw=lw, solid_capstyle="butt")
ax.text(XR, 0.932, "WHO SEES IT", fontsize=10.5, color=GREY, va="center")
rules(0.910, RULE, 1.0)
rows = [
 (0.872, "1", "Outcome observed",
  [r"Endpoint information $Y_i$ arrives from the trial sites;", r"the assignment $T_i$ stays sealed in the", "randomization system."],
  ["trial sites", "monitoring system"]),
 (0.680, "2", "Wager recorded",
  [r"$\lambda_i$ is computed from $Y_i$ and the history $\mathcal{F}_{i-1}$, then", "recorded before the assignment is revealed."],
  ["monitoring system"]),
 (0.527, "3", "Assignment revealed",
  [r"$T_i$ is then released, solely to the independent", "unblinded monitoring system."],
  ["monitoring system only"]),
 (0.374, "4", "Wealth updated",
  [r"$W_i = W_{i-1} \times \lambda_i / p_i$ if $T_i = 1$, else $\times\,(1-\lambda_i)/(1-p_i)$;", "compared with the threshold."],
  ["monitoring system"]),
]
LH = 0.0385
for y, num, head, body, who in rows:
    ax.text(XN, y - 0.012, num, fontsize=30, color="#3a3a3a", va="center")
    ax.text(XT, y, head, fontsize=13.5, color=INK, weight="bold", va="center")
    for k, line in enumerate(body): ax.text(XT, y - 0.036 - k * LH, line, fontsize=12.5, color=INK, va="center")
    for k, w in enumerate(who): ax.text(XR, y - k * 0.031, w, fontsize=10.8, color=GREY, style="italic", va="center")
for y in (0.718, 0.565, 0.412): rules(y, SOFT, 1.0)
# wealth strip
x0, y0, w, h = XN, 0.045, XE - XN, 0.195
ax.add_patch(Rectangle((x0, y0), w, h, fc="#efeeea", ec="#dcdcdc", lw=1.0))
thr = y0 + 0.148
ax.plot([x0 + 0.03, x0 + w - 0.035], [thr, thr], color="#444444", lw=1.1, dashes=(4, 4))
ys = [0.075, 0.078, 0.073, 0.084, 0.081, 0.093, 0.090, 0.103, 0.099, 0.114, 0.110, 0.129, 0.138, 0.153]
xs = [x0 + 0.03 + k * (w - 0.09) / (len(ys) - 1) for k in range(len(ys))]
ax.plot(xs, [y0 + v for v in ys], color=INK, lw=1.6, solid_joinstyle="round")
ax.plot(xs[-1], y0 + ys[-1], "o", color=INK, ms=4.5)
ax.text(x0 + 0.03, y0 + 0.030, r"wealth $W_i$, inspectable at any update", fontsize=10.5, color=GREY, style="italic", va="center")
ax.text(x0 + w - 0.03, y0 + 0.030, r"threshold $1/\alpha = 20$", fontsize=10.5, color=GREY, style="italic", va="center", ha="right")
ax.plot([XN, XE], [0.022, 0.022], color=RULE, lw=1.0)
for ext in ("png", "pdf"): fig.savefig(out / f"fig1_infoflow.{ext}", dpi=300)
print("written:", out / "fig1_infoflow.png")
