import matplotlib.pyplot as plt
import matplotlib.ticker as mticker
import scienceplots

# Use SciencePlots style; add 'no-latex' if you do not have LaTeX installed
plt.style.use(["science", "grid"])  # or ["science", "grid", "no-latex"]

# -----------------------
# Data
# -----------------------
batch = [1, 8, 32, 128, 512, 2048]
qps_cube = [9.225e5, 9.291e5, 9.272e5, 9.324e5, 9.323e5, 9.320e5]
qps_cyl  = [2.546e5, 2.547e5, 2.552e5, 2.520e5, 2.553e5, 2.551e5]

defl = [1e-2, 5e-3, 1e-3, 5e-4]
tess_time_cyl = [8.88e-4, 5.34e-4, 7.899e-3, 1.2694e-2]
tris_cyl = [1084, 2108, 16636, 33020]

# -----------------------
# Helper: consistent save
# -----------------------
def savefig(fig, stem, dpi=600):
    fig.tight_layout()
    fig.savefig(f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(f"{stem}.png", dpi=dpi, bbox_inches="tight")

# -----------------------
# 1) Triangle count vs deflection (cylinder)
# -----------------------
fig, ax = plt.subplots(figsize=(3.4, 2.6))  # ~single-column figure
ax.set_xscale("log")
ax.set_yscale("log")

ax.plot(defl, tris_cyl, marker="o", linewidth=1.6, markersize=4)

ax.set_xlabel("Chordal deflection")
ax.set_ylabel("Triangle count")
ax.set_title("Cylinder tessellation density")

# Log tick formatting: show 10^n
ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
ax.xaxis.set_minor_locator(mticker.LogLocator(base=10, subs="auto"))
ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs="auto"))

savefig(fig, "fig_triangle_count_vs_deflection")

# -----------------------
# 2) Tessellation time vs deflection (cylinder)
# -----------------------
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.set_xscale("log")
ax.set_yscale("log")

ax.plot(defl, tess_time_cyl, marker="o", linewidth=1.6, markersize=4)

ax.set_xlabel("Chordal deflection")
ax.set_ylabel("Tessellation time (s)")
ax.set_title("Cylinder tessellation cost")

ax.xaxis.set_major_locator(mticker.LogLocator(base=10))
ax.xaxis.set_minor_locator(mticker.LogLocator(base=10, subs="auto"))
ax.yaxis.set_major_locator(mticker.LogLocator(base=10))
ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs="auto"))

savefig(fig, "fig_tessellation_time_vs_deflection")

# -----------------------
# 3) Inverse-eval throughput vs batch size (cube vs cylinder)
# -----------------------
fig, ax = plt.subplots(figsize=(3.4, 2.6))
ax.set_xscale("log")

ax.plot(batch, qps_cube, marker="o", linewidth=1.6, markersize=4, label="Cube")
ax.plot(batch, qps_cyl,  marker="o", linewidth=1.6, markersize=4, label="Cylinder")

ax.set_xlabel("Batch size")
ax.set_ylabel(r"Inverse-eval throughput (q s$^{-1}$)")
ax.set_title("Inverse evaluation throughput")

# Make y-axis labels scientific and clean
ax.yaxis.set_major_formatter(mticker.ScalarFormatter(useMathText=True))
ax.ticklabel_format(axis="y", style="sci", scilimits=(0, 0))

ax.legend(frameon=True, fontsize=8)

savefig(fig, "fig_inveval_throughput_vs_batch")

plt.show()
