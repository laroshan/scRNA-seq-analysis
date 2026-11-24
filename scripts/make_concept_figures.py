"""Draw the teaching figures used in CONCEPTS.md.

Every figure uses small, made-up toy data, so it explains an idea and is not a result.
Run from the project root:

    python scripts/make_concept_figures.py

Writes PNG files to figures/concepts/ and the toy-example numbers to figures/concepts/numbers.md.
Needs only numpy, scipy and matplotlib. Written to run on Python 3.8 or newer.
"""

import os
from itertools import combinations

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402
from scipy.optimize import linear_sum_assignment  # noqa: E402
from scipy.spatial import cKDTree  # noqa: E402

OUT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "figures", "concepts")

# Palette: the first three categorical slots (validated all-pairs, light surface) plus text and surface tokens.
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, MUTED, GRID, SURFACE = "#0b0b0b", "#52514e", "#8a8984", "#e4e3df", "#fcfcfb"
SEQ = LinearSegmentedColormap.from_list("seq_blue", ["#ffffff", "#cde2fb", "#86b6ef", "#3987e5", "#1c5cab", "#0d366b"])
DIV = LinearSegmentedColormap.from_list("div", ["#1c5cab", "#86b6ef", "#f0efec", "#f0a3a2", "#c03636"])

plt.rcParams.update({
    "figure.facecolor": SURFACE, "axes.facecolor": SURFACE, "savefig.facecolor": SURFACE,
    "font.size": 11, "axes.titlesize": 12, "axes.titleweight": "bold", "axes.titlelocation": "left",
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8, "axes.axisbelow": True,
    "legend.frameon": False, "lines.linewidth": 2,
})


def save(fig, name):
    fig.savefig(os.path.join(OUT, name), dpi=150, bbox_inches="tight")
    plt.close(fig)
    print("wrote", name)


def blobs(rng, centres, sizes, sd):
    pts, lab = [], []
    for i, (c, n) in enumerate(zip(centres, sizes)):
        sdi = sd[i] if isinstance(sd, (list, tuple)) else sd
        pts.append(rng.normal(c, sdi, size=(n, 2)))
        lab += [i] * n
    return np.vstack(pts), np.array(lab)


# ---------- metrics written out by hand, so the reader can follow them ----------
def comb2(x):
    return x * (x - 1) / 2.0


def contingency(a, b):
    ua, ub = np.unique(a), np.unique(b)
    return np.array([[np.sum((a == i) & (b == j)) for j in ub] for i in ua], dtype=float)


def ari(a, b):
    c = contingency(a, b)
    n = c.sum()
    s_ij = comb2(c).sum()
    s_a, s_b = comb2(c.sum(1)).sum(), comb2(c.sum(0)).sum()
    expected = s_a * s_b / comb2(n)
    return (s_ij - expected) / (0.5 * (s_a + s_b) - expected)


def nmi(a, b):
    c = contingency(a, b)
    p = c / c.sum()
    pa, pb = p.sum(1), p.sum(0)
    nz = p > 0
    mi = np.sum(p[nz] * np.log(p[nz] / np.outer(pa, pb)[nz]))
    h = lambda q: -np.sum(q[q > 0] * np.log(q[q > 0]))  # noqa: E731
    return mi / (0.5 * (h(pa) + h(pb)))


def hungarian_accuracy(a, b):
    c = contingency(a, b)
    r, k = linear_sum_assignment(-c)
    return c[r, k].sum() / c.sum()


def best_match_f1(a, b, label):
    t = a == label
    best = 0.0
    for cl in np.unique(b):
        c = b == cl
        both = np.sum(t & c)
        if both:
            best = max(best, 2 * both / (t.sum() + c.sum()))
    return best


# ---------- figures ----------
def fig_cells_as_points(rng):
    X, y = blobs(rng, [(6, 0.6), (0.6, 5.5), (0.8, 0.8)], [25, 25, 25], 0.55)
    X = np.clip(X, 0, None)
    fig, ax = plt.subplots(figsize=(6.4, 4.6))
    names = ["T cells", "B cells", "Monocytes"]
    for i, col in enumerate([BLUE, ORANGE, AQUA]):
        ax.scatter(X[y == i, 0], X[y == i, 1], s=38, color=col, edgecolor=SURFACE, linewidth=1, label=names[i])
    ax.annotate("one cell = one point", xy=X[3], xytext=(3.2, 3.2), color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=1))
    ax.set_xlabel("expression of CD3E (a T-cell gene)")
    ax.set_ylabel("expression of MS4A1 (a B-cell gene)")
    ax.set_title("Each cell is a point; similar cells sit close together")
    ax.legend(loc="upper right")
    save(fig, "01_cells_as_points.png")


def fig_distance():
    a, b = np.array([1.0, 4.0]), np.array([5.0, 1.0])
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.2))
    ax = axes[0]
    ax.plot([a[0], b[0]], [a[1], b[1]], color=INK2, lw=1.5, ls="-")
    ax.plot([a[0], b[0], b[0]], [a[1], a[1], b[1]], color=MUTED, lw=1, ls=":")
    for p, n, col in [(a, "cell A", BLUE), (b, "cell B", ORANGE)]:
        ax.scatter(*p, s=80, color=col, zorder=3)
        ax.annotate(n, p, xytext=(8, 6), textcoords="offset points")
    ax.text(0.3, 1.0, "distance = √(4² + 3²) = 5", color=INK)
    ax.text(2.6, 4.15, "4 apart on gene 1", color=INK2, fontsize=10)
    ax.text(5.1, 2.4, "3 apart\non gene 2", color=INK2, fontsize=10)
    ax.set_xlim(0, 7); ax.set_ylim(0, 5.5)
    ax.set_xlabel("gene 1"); ax.set_ylabel("gene 2")
    ax.set_title("Straight-line (Euclidean) distance")

    ax = axes[1]
    c, d = np.array([1.0, 2.0]), np.array([3.0, 6.0])
    e = np.array([5.0, 1.0])
    for p, n, col in [(c, "cell C", BLUE), (d, "cell D (C × 3)", BLUE), (e, "cell E", ORANGE)]:
        ax.annotate("", xy=p, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=col, lw=2))
        ax.annotate(n, p, xytext=(6, 4), textcoords="offset points")
    ax.set_xlim(0, 7); ax.set_ylim(0, 7)
    ax.set_xlabel("gene 1"); ax.set_ylabel("gene 2")
    ax.set_title("Cosine: compares direction, ignores length")
    ax.text(3.4, 3.3, "C and D point the same way:\ncosine distance 0", color=INK2, fontsize=10)
    save(fig, "02_distance.png")


def fig_spread(rng):
    g1 = rng.normal(5, 0.3, 40)
    g2 = rng.normal(5, 1.6, 40)
    fig, ax = plt.subplots(figsize=(8, 2.8))
    for row, vals, name, col in [(1, g1, "gene X: almost the same in every cell", MUTED),
                                 (0, g2, "gene Y: varies a lot between cells", BLUE)]:
        ax.scatter(vals, np.full_like(vals, row) + rng.uniform(-0.08, 0.08, vals.size), s=30, color=col,
                   edgecolor=SURFACE, linewidth=0.8)
        ax.plot([vals.mean()] * 2, [row - 0.25, row + 0.25], color=INK, lw=2)
        ax.text(0.2, row + 0.3, f"{name}   (variance {vals.var():.2f})", color=INK2, fontsize=10)
    ax.set_yticks([]); ax.set_ylim(-0.5, 1.7); ax.set_xlim(0, 10)
    ax.set_xlabel("expression level (black tick = mean)")
    ax.set_title("Variance = how spread out the values are")
    ax.grid(axis="y", visible=False)
    save(fig, "03_variance.png")


def fig_regression_vs_pca(rng):
    x = rng.normal(0, 1.6, 40)
    y = 0.7 * x + rng.normal(0, 0.7, 40)
    X = np.c_[x, y] - np.c_[x, y].mean(0)
    slope, icpt = np.polyfit(X[:, 0], X[:, 1], 1)
    _, _, vt = np.linalg.svd(X, full_matrices=False)
    v = vt[0]
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4), sharex=True, sharey=True)
    xs = np.linspace(-4.5, 4.5, 2)
    ax = axes[0]
    for px, py in X:
        ax.plot([px, px], [py, slope * px + icpt], color=ORANGE, lw=1)
    ax.plot(xs, slope * xs + icpt, color=INK, lw=2)
    ax.scatter(X[:, 0], X[:, 1], s=30, color=BLUE, zorder=3, edgecolor=SURFACE, linewidth=0.8)
    ax.set_title("Regression: predict gene 2 from gene 1")
    ax.text(-4.4, 3.3, "minimises the VERTICAL gaps (orange)\ngene 2 is the 'answer'", color=INK2, fontsize=10)
    ax = axes[1]
    for p in X:
        q = (p @ v) * v
        ax.plot([p[0], q[0]], [p[1], q[1]], color=ORANGE, lw=1)
    ax.plot(xs * v[0] / abs(v[0]) * 1.2, xs * v[1] / abs(v[0]) * 1.2, color=INK, lw=2)
    ax.scatter(X[:, 0], X[:, 1], s=30, color=BLUE, zorder=3, edgecolor=SURFACE, linewidth=0.8)
    ax.set_title("PCA: the line through the spread")
    ax.text(-4.4, 3.3, "minimises the PERPENDICULAR gaps\nno gene is the 'answer'", color=INK2, fontsize=10)
    for ax in axes:
        ax.set_xlim(-4.6, 4.6); ax.set_ylim(-4, 4.2); ax.set_aspect("equal")
        ax.set_xlabel("gene 1 (centred)")
    axes[0].set_ylabel("gene 2 (centred)")
    save(fig, "04_regression_vs_pca.png")


def fig_pca_steps(rng):
    t = rng.normal(0, 1.8, 60)
    X = np.c_[t + 6, 0.55 * t + 4 + rng.normal(0, 0.45, 60)]
    m = X.mean(0)
    Xc = X - m
    _, s, vt = np.linalg.svd(Xc, full_matrices=False)
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.6), gridspec_kw={"width_ratios": [1, 1.2, 1], "wspace": 0.3})
    ax = axes[0]
    ax.scatter(X[:, 0], X[:, 1], s=26, color=BLUE, edgecolor=SURFACE, linewidth=0.7)
    ax.scatter(*m, s=110, marker="X", color=INK, zorder=4)
    ax.annotate("mean of each gene", m, xytext=(m[0] + 0.6, m[1] - 2.2), color=INK2,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=1))
    ax.set_xlim(0, 12); ax.set_ylim(0, 8)
    ax.set_title("1. Raw data"); ax.set_xlabel("gene 1"); ax.set_ylabel("gene 2")
    ax = axes[1]
    ax.scatter(Xc[:, 0], Xc[:, 1], s=26, color=BLUE, edgecolor=SURFACE, linewidth=0.7)
    sign = np.sign(vt[:, 0:1]) * np.array([[1.0], [-1.0]])  # PC1 points right, PC2 points up-left
    for i, (lab, col, off) in enumerate([("PC1: most spread", ORANGE, (-3.2, 0.6)), ("PC2: the rest", AQUA, (-2.6, 0.5))]):
        vec = vt[i] * sign[i] * s[i] / np.sqrt(len(X)) * 2.5
        ax.annotate("", xy=vec, xytext=(0, 0), arrowprops=dict(arrowstyle="->", color=col, lw=2.5))
        ax.text(vec[0] + off[0], vec[1] + off[1], lab, color=INK, fontsize=10)
    ax.set_xlim(-6, 6); ax.set_ylim(-4, 4); ax.set_aspect("equal")
    ax.set_title("2. Centre, then find the directions"); ax.set_xlabel("gene 1 − mean"); ax.set_ylabel("gene 2 − mean")
    ax = axes[2]
    pc1 = Xc @ vt[0]
    ax.scatter(pc1, np.zeros_like(pc1) + rng.uniform(-0.05, 0.05, pc1.size), s=26, color=BLUE,
               edgecolor=SURFACE, linewidth=0.7)
    ax.set_ylim(-1, 1); ax.set_yticks([])
    ax.set_title("3. Keep PC1 only: 2 numbers → 1"); ax.set_xlabel("position along PC1")
    ax.grid(axis="y", visible=False)
    share = s[0] ** 2 / np.sum(s ** 2)
    ax.text(pc1.min(), 0.45, f"PC1 keeps {share:.0%} of the spread", color=INK2, fontsize=10)
    save(fig, "05_pca_steps.png")


def toy_counts(rng, n_cells=30, n_genes=40):
    lam = np.full((n_cells, n_genes), 0.05)
    groups = np.repeat([0, 1, 2], n_cells // 3)
    for g in range(3):
        lam[groups == g, g * 10:(g * 10) + 10] = 2.5
    return rng.poisson(lam).astype(float), groups


def fig_sparse_centring(rng):
    C, _ = toy_counts(rng)
    L = np.log1p(C)
    Lc = L - L.mean(0)
    fig, axes = plt.subplots(1, 2, figsize=(13, 4.6), gridspec_kw={"wspace": 0.35})
    ax = axes[0]
    im = ax.imshow(L, cmap=SEQ, aspect="auto", vmin=0, interpolation="nearest")
    ax.set_title(f"Log counts: {np.mean(L == 0):.0%} of entries are exactly 0 (white)")
    ax.set_xlabel("genes"); ax.set_ylabel("cells"); ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.04, label="log(1 + count)")
    ax = axes[1]
    lim = np.abs(Lc).max()
    im = ax.imshow(Lc, cmap=DIV, aspect="auto", vmin=-lim, vmax=lim, interpolation="nearest")
    ax.set_title(f"After subtracting each gene's mean: {np.mean(np.abs(Lc) < 1e-12):.0%} are 0")
    ax.set_xlabel("genes"); ax.grid(False)
    fig.colorbar(im, ax=ax, fraction=0.04, label="value − gene mean")
    save(fig, "06_sparse_then_centred.png")


def fig_memory():
    n = np.array([10_000, 25_000, 50_000, 100_000, 250_000], dtype=float)
    genes = 2000
    sparse = n * genes * 0.10 * 8 / 1e9          # 10% non-zero, 4-byte value + 4-byte column index
    dense = n * genes * 4 / 1e9                  # float32
    pdist = n * (n - 1) / 2 * 8 / 1e9            # every pair, float64
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    for vals, lab, col in [(pdist, "all pairwise distances (hierarchical, Dunn)", ORANGE),
                           (dense, "dense copy (needed by full / randomized PCA)", BLUE),
                           (sparse, "sparse matrix (assuming 10% non-zero)", AQUA)]:
        ax.plot(n, vals, marker="o", markersize=7, color=col, markeredgecolor=SURFACE, markeredgewidth=1.5)
        ax.text(n[-1] * 1.08, vals[-1], f"{lab}\n{vals[-1]:,.1f} GB at 250k", va="center", fontsize=10, color=INK)
    ax.axhline(24, color=INK2, lw=1)
    ax.text(n[0], 30, "your Mac: 24 GB", color=INK2, fontsize=10)
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xticks(n); ax.set_xticklabels(["10k", "25k", "50k", "100k", "250k"])
    ax.minorticks_off()
    ax.set_xlabel("number of cells (log scale)"); ax.set_ylabel("memory, GB (log scale)")
    ax.set_title("Memory grows with cells, or with cells² for pairwise distances")
    save(fig, "07_memory_growth.png")


def toy_highdim(rng, n=600, d=2000):
    lam = np.full((n, d), 0.08)
    groups = rng.integers(0, 3, n)
    for g in range(3):
        genes = rng.choice(d, 120, replace=False)
        lam[np.ix_(groups == g, genes)] = 1.5
    return np.log1p(rng.poisson(lam)).astype(float), groups


def fig_rp_ratios(rng):
    X, _ = toy_highdim(rng)
    pairs = rng.integers(0, X.shape[0], size=(3000, 2))
    pairs = pairs[pairs[:, 0] != pairs[:, 1]]
    orig = np.linalg.norm(X[pairs[:, 0]] - X[pairs[:, 1]], axis=1)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8), sharey=True)
    for ax, k in zip(axes, [10, 50, 500]):
        R = rng.normal(0, 1 / np.sqrt(k), size=(X.shape[1], k))
        P = X @ R
        ratio = np.linalg.norm(P[pairs[:, 0]] - P[pairs[:, 1]], axis=1) / orig
        ax.hist(ratio, bins=np.linspace(0.3, 1.7, 57), color=BLUE, edgecolor=SURFACE, linewidth=0.6)
        ax.axvline(1, color=INK, lw=1.5)
        lo, hi = np.percentile(ratio, [5, 95])
        ax.set_title(f"→ {k} dimensions")
        ax.set_xlabel("new distance ÷ old distance")
        ax.text(0.04, 0.92, f"90% of pairs: {lo:.2f}–{hi:.2f}", transform=ax.transAxes, color=INK2, fontsize=10)
    axes[0].set_ylabel("number of cell pairs")
    axes[0].set_ylim(0, 1150)
    fig.suptitle("Random projection of 600 toy cells × 2,000 genes: how much each distance changed",
                 x=0.01, ha="left", fontweight="bold", y=1.04)
    save(fig, "08_rp_keeps_distances.png")


def fig_rp_matrices(rng):
    k, d = 16, 60
    X, _ = toy_counts(rng, 30, d)
    data_density = np.mean(X != 0)
    gauss = rng.normal(0, 1, size=(k, d))
    s = np.sqrt(d)
    sparse = rng.choice([-1.0, 0.0, 1.0], size=(k, d), p=[1 / (2 * s), 1 - 1 / s, 1 / (2 * s)])
    match = np.where(rng.random((k, d)) < data_density, rng.normal(0, 1, (k, d)), 0.0)
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.7))
    for ax, M, title in [(axes[0], gauss, "Gaussian RP\nevery entry random"),
                         (axes[1], sparse, f"Sparse RP\ndensity 1/√{d} ≈ {1 / s:.0%}"),
                         (axes[2], match, f"Matching Sparsity RP\ndensity = the data's ({data_density:.0%})")]:
        lim = np.abs(M).max()
        ax.imshow(M, cmap=DIV, vmin=-lim, vmax=lim, aspect="auto", interpolation="nearest")
        ax.set_title(title, fontsize=11)
        ax.set_xlabel(f"genes (input) · {np.mean(M != 0):.0%} of entries non-zero"); ax.set_yticks([]); ax.grid(False)
    axes[0].set_ylabel("new dimensions (output)")
    save(fig, "09_three_random_matrices.png")


def fig_rare_population(rng):
    X, y = blobs(rng, [(-3, 0), (3, 0), (3.3, 1.9)], [200, 200, 8], [(1.1, 0.45), (1.1, 0.45), (0.25, 0.2)])
    Xc = X - X.mean(0)
    _, _, vt = np.linalg.svd(Xc, full_matrices=False)
    pc1 = Xc @ vt[0]
    names = ["common type A", "common type B", "rare type (2%)"]
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.2), gridspec_kw={"width_ratios": [1.2, 1]})
    ax = axes[0]
    for i, col in enumerate([BLUE, ORANGE, AQUA]):
        ax.scatter(X[y == i, 0], X[y == i, 1], s=30 if i < 2 else 55, color=col, edgecolor=SURFACE, linewidth=0.8,
                   zorder=3 if i == 2 else 2)
    ax.annotate("", xy=vt[0] * 5 + X.mean(0), xytext=-vt[0] * 5 + X.mean(0),
                arrowprops=dict(arrowstyle="<->", color=INK, lw=1.5))
    ax.text(-6.2, 1.1, "PC1 (most spread)", color=INK, fontsize=10)
    ax.annotate(names[2], X[y == 2].mean(0), xytext=(0.0, 2.9), color=INK,
                arrowprops=dict(arrowstyle="->", color=INK2, lw=1))
    ax.text(-6.3, -2.1, "common type A", color=INK2, fontsize=10)
    ax.text(2.2, -2.1, "common type B", color=INK2, fontsize=10)
    ax.set_xlim(-7, 7); ax.set_ylim(-2.6, 3.6)
    ax.set_title("2 genes: the rare type is clearly separate")
    ax.set_xlabel("gene 1"); ax.set_ylabel("gene 2")
    ax = axes[1]
    for i, col in enumerate([BLUE, ORANGE, AQUA]):
        ax.scatter(pc1[y == i], np.full((y == i).sum(), 0.0) + (0.25 if i == 2 else 0) + rng.uniform(-0.06, 0.06, (y == i).sum()),
                   s=30 if i < 2 else 55, color=col, edgecolor=SURFACE, linewidth=0.8)
    ax.text(pc1[y == 2].mean() - 1.0, 0.45, "rare type lands\ninside type B", color=INK, fontsize=10)
    ax.set_ylim(-0.6, 1.0); ax.set_yticks([]); ax.grid(axis="y", visible=False)
    ax.set_title("Kept only PC1: the rare type is hidden")
    ax.set_xlabel("position along PC1 (rare cells drawn higher to be seen)")
    save(fig, "10_rare_type_and_pc1.png")


def fig_kmeans(rng):
    X, _ = blobs(rng, [(1, 1), (5, 1.5), (3, 5)], [30, 30, 30], 0.7)
    C = np.array([[0.5, 4.5], [1.5, 0.5], [2.5, 1.0]])
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    cols = np.array([BLUE, ORANGE, AQUA])
    titles = ["1. Drop k = 3 random centres", "2. Each cell joins its nearest centre", "3. Move centres, repeat until stable"]
    for step, ax in enumerate(axes):
        if step == 0:
            ax.scatter(X[:, 0], X[:, 1], s=26, color=MUTED, edgecolor=SURFACE, linewidth=0.7)
            centres = C
        else:
            if step == 2:
                for _ in range(10):
                    lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), 1)
                    C = np.array([X[lab == j].mean(0) for j in range(3)])
            lab = np.argmin(((X[:, None, :] - C[None]) ** 2).sum(-1), 1)
            ax.scatter(X[:, 0], X[:, 1], s=26, c=cols[lab], edgecolor=SURFACE, linewidth=0.7)
            centres = C
        ax.scatter(centres[:, 0], centres[:, 1], s=180, marker="X", c=list(cols), edgecolor=INK, linewidth=1.2, zorder=4)
        ax.set_title(titles[step], fontsize=11); ax.set_xlim(-1, 7); ax.set_ylim(-1, 7)
    save(fig, "11_kmeans_steps.png")


def fig_spherical(rng):
    X, y = blobs(rng, [(1.5, 0.4), (4.5, 1.2), (0.5, 2.5)], [12, 12, 12], 0.3)
    U = X / np.linalg.norm(X, axis=1, keepdims=True)
    fig, axes = plt.subplots(1, 2, figsize=(10, 4.4))
    for ax, P, title in [(axes[0], X, "Before: points at different distances"),
                         (axes[1], U, "After scaling each cell to length 1")]:
        for i, col in enumerate([BLUE, ORANGE, AQUA]):
            ax.scatter(P[y == i, 0], P[y == i, 1], s=30, color=col, edgecolor=SURFACE, linewidth=0.8)
        ax.set_aspect("equal"); ax.set_title(title, fontsize=11)
    th = np.linspace(0, np.pi / 2, 100)
    axes[1].plot(np.cos(th), np.sin(th), color=MUTED, lw=1)
    axes[1].set_xlim(-0.1, 1.2); axes[1].set_ylim(-0.3, 1.2)
    axes[1].text(0.05, 1.08, "blue and orange now overlap: same direction", color=INK2, fontsize=10)
    axes[0].set_xlim(-0.5, 6); axes[0].set_ylim(-0.5, 4)
    save(fig, "12_spherical_kmeans.png")


def fig_knn_graph(rng):
    X, y = blobs(rng, [(1, 1), (5, 1.3), (3, 4.6)], [22, 22, 22], 0.75)
    tree = cKDTree(X)
    _, idx = tree.query(X, k=5)
    fig, ax = plt.subplots(figsize=(6.6, 5))
    for i in range(len(X)):
        for j in idx[i, 1:]:
            ax.plot([X[i, 0], X[j, 0]], [X[i, 1], X[j, 1]], color=MUTED, lw=0.7, zorder=1)
    for i, col in enumerate([BLUE, ORANGE, AQUA]):
        ax.scatter(X[y == i, 0], X[y == i, 1], s=40, color=col, edgecolor=SURFACE, linewidth=1, zorder=3)
    ax.set_title("Neighbour graph: each cell linked to its 4 nearest cells")
    ax.text(-0.9, 6.3, "Leiden looks for groups with many links inside, few between", color=INK2, fontsize=10)
    ax.set_xlim(-1, 7.5); ax.set_ylim(-1, 6.8); ax.grid(False)
    save(fig, "13_knn_graph.png")


def fig_ari_hides_rare():
    truth = np.array([0] * 600 + [1] * 380 + [2] * 20)
    clusters = truth.copy()
    clusters[truth == 2] = 1          # the rare type is merged into the big type B
    clusters[:15] = 1                 # and a few A cells are wrong too
    scores = [("ARI", ari(truth, clusters)), ("NMI", nmi(truth, clusters)),
              ("Hungarian accuracy", hungarian_accuracy(truth, clusters)),
              ("F1, rare type", best_match_f1(truth, clusters, 2))]
    fig, ax = plt.subplots(figsize=(7.6, 3.4))
    names = [s[0] for s in scores][::-1]
    vals = [s[1] for s in scores][::-1]
    cols = [ORANGE] + [BLUE] * 3
    ax.barh(names, vals, color=cols, height=0.55)
    for i, v in enumerate(vals):
        ax.text(v + 0.02, i, f"{v:.2f}", va="center", color=INK)
    ax.set_xlim(0, 1.12); ax.set_xlabel("score (1 = perfect)")
    ax.set_title("1,000 cells; the 20 rare cells were merged into a big group")
    ax.grid(axis="y", visible=False)
    save(fig, "14_ari_hides_rare.png")
    return scores


def fig_seeds(rng):
    fig, axes = plt.subplots(1, 2, figsize=(10.5, 3.6), sharey=True)
    cases = [("Ranges overlap → 'no clear difference'", [0.69, 0.70, 0.715, 0.72, 0.735], [0.675, 0.69, 0.70, 0.71, 0.725]),
             ("Ranges apart → a real difference", [0.765, 0.77, 0.775, 0.78, 0.79], [0.645, 0.65, 0.66, 0.665, 0.675])]
    for ax, (title, a, b) in zip(axes, cases):
        for x, vals, col in [(0, np.array(a), BLUE), (1, np.array(b), ORANGE)]:
            ax.scatter(np.full(5, x) + rng.uniform(-0.06, 0.06, 5), vals, s=50, color=col, edgecolor=SURFACE, zorder=3)
            ax.plot([x - 0.2, x + 0.2], [np.median(vals)] * 2, color=INK, lw=2)
            ax.plot([x, x], [vals.min(), vals.max()], color=INK2, lw=1)
        ax.set_xticks([0, 1]); ax.set_xticklabels(["method 1", "method 2"])
        ax.set_xlim(-0.6, 1.6); ax.set_title(title, fontsize=11); ax.grid(axis="x", visible=False)
    axes[0].set_ylabel("ARI (made-up numbers)")
    axes[0].set_ylim(0.6, 0.85)
    fig.text(0.01, -0.04, "Dots = 5 seeds. Black bar = median. Grey line = min–max range.", color=INK2, fontsize=10)
    save(fig, "15_seeds_and_spread.png")


def write_numbers(scores):
    labels = np.array(list("TTTTBBBMMM"))
    clusters = np.array([1, 1, 1, 2, 2, 2, 2, 3, 3, 3])
    lines = ["# Numbers used in CONCEPTS.md (made-up toy examples)", "",
             "## 10-cell example", "",
             "labels:   " + " ".join(labels), "clusters: " + " ".join(map(str, clusters)), "",
             f"- ARI = {ari(labels, clusters):.3f}",
             f"- NMI = {nmi(labels, clusters):.3f}",
             f"- Hungarian accuracy = {hungarian_accuracy(labels, clusters):.3f}"]
    for t in "TBM":
        lines.append(f"- best-match F1, {t} = {best_match_f1(labels, clusters, t):.3f}")
    same_lab = sum(1 for i, j in combinations(range(10), 2) if labels[i] == labels[j])
    same_cl = sum(1 for i, j in combinations(range(10), 2) if clusters[i] == clusters[j])
    both = sum(1 for i, j in combinations(range(10), 2) if labels[i] == labels[j] and clusters[i] == clusters[j])
    lines += ["", f"- pairs of cells: 45; same label: {same_lab}; same cluster: {same_cl}; same in both: {both}"]
    lines += ["", "## Rare type merged (figure 14)", ""] + [f"- {n} = {v:.3f}" for n, v in scores]
    with open(os.path.join(OUT, "numbers.md"), "w") as f:
        f.write("\n".join(lines) + "\n")
    print("\n".join(lines))


def main():
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(2026)
    fig_cells_as_points(rng)
    fig_distance()
    fig_spread(rng)
    fig_regression_vs_pca(rng)
    fig_pca_steps(rng)
    fig_sparse_centring(rng)
    fig_memory()
    fig_rp_ratios(rng)
    fig_rp_matrices(rng)
    fig_rare_population(rng)
    fig_kmeans(rng)
    fig_spherical(rng)
    fig_knn_graph(rng)
    scores = fig_ari_hides_rare()
    fig_seeds(rng)
    write_numbers(scores)


if __name__ == "__main__":
    main()
