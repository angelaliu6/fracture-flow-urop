"""Generate the Section 4 figure set comparing the flow-model surrogate
("Updated Code/flow distrib.ipynb") and the transport-model surrogate
("Updated Code/[per p, t] probability of >=p percent having bkt <= t.ipynb").
Run standalone: `.venv/bin/python "Updated Code/figures/generate_figures.py"`.

The transport-model pipeline (model_constrained, error_score_df, val_ps,
and the plotting helpers) isn't saved as a loadable artifact anywhere --
only the .keras weights are -- so it's regenerated here by re-running the
relevant prefix of that notebook's own cells in-process (same values each
time: no randomness in the data/normalization, and retraining the NNs
takes ~15-20s and gives an equivalent-quality model on every run). This
keeps a single source of truth (the notebook) instead of a second,
driftable copy of that logic living in this script.

Figure set, matching the restructured Section 4 (flow story, then transport):
    fig1_aggregate_accuracy.png   held-out error histograms, both models
    fig2a_parity_flow.png         mixture-coefficient parity (flow)
    fig2b_parity_transport.png    pooled F parity (transport)
    fig6_error_vs_parameters.png  held-out error vs each theta component
    fig7_sweeps_flow.png          flow response to r, d, b  (1x3)
    fig8_sweeps_transport.png     transport response to r, d, b  (1x3)
    figA1_lambda_sweep.png        constraint-weight sweep (appendix)

NOT generated here: fig0_flow_examples.png comes from flow distrib.ipynb.
Adding the directly-fitted GMM alongside the surrogate's prediction (so the
figure can show family adequacy separately from surrogate error) has to be
done in that notebook.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.stats import norm as scipy_norm

# --- Page-aware type scaling -------------------------------------------------
# These figures are rendered several inches wider than the column they are
# printed in, so LaTeX shrinks them and the type shrinks with it: at the old
# 20 in render width, a 19 pt label arrived on the page at about 5.7 pt.
# Raising rcParams alone cannot fix this -- the render size is the lever.
# So: the figsizes below are much closer to the printed size, and TARGET_PT
# specifies the sizes we want AFTER the remaining shrink. scale_fonts_for()
# multiplies TARGET_PT by each figure's own downscale factor, so every figure
# lands at the same on-page type size regardless of how wide it was rendered.
#
# Verify TEXTWIDTH_IN by putting \the\textwidth in the manuscript and
# dividing the printed value (in pt) by 72.27.
TEXTWIDTH_IN = 6.0

TARGET_PT = {                        # desired ON-PAGE point sizes
    "font.size":              10.0,  # manuscript body text is 10 pt
    "axes.labelsize":         10.0,
    "axes.titlesize":         10.0,
    "xtick.labelsize":         9.0,
    "ytick.labelsize":         9.0,
    "legend.fontsize":         9.0,
    "legend.title_fontsize":   9.0,
    "figure.titlesize":       11.5,
}


def scale_fonts_for(fig_width_in, frac=1.0):
    """Set rcParams so text renders at TARGET_PT once a figure of width
    `fig_width_in` inches is included at `frac` * \textwidth. Returns the
    scale factor, for sizing scatter markers (whose `s` is in pt^2)."""
    k = fig_width_in / (frac * TEXTWIDTH_IN)
    rc = {key: pt * k for key, pt in TARGET_PT.items()}
    rc.update({                      # keep strokes proportionate to the type
        "axes.linewidth":      0.8 * k,
        "xtick.major.width":   0.8 * k,
        "ytick.major.width":   0.8 * k,
        "xtick.major.size":    3.5 * k,
        "ytick.major.size":    3.5 * k,
        "lines.linewidth":     1.5 * k,
        "lines.markersize":    6.0 * k,
        "patch.linewidth":     1.0 * k,
        "legend.frameon":      False,
        "figure.dpi":          150,
    })
    matplotlib.rcParams.update(rc)
    return k


FIGURES_DIR = Path(__file__).parent
UPDATED_CODE_DIR = FIGURES_DIR.parent

# Color convention shared by every figure: black/gray = ground-truth /
# empirical / simulated data, tab:blue = flow-model NN surrogate,
# tab:orange = transport-model NN surrogate, viridis = continuous
# parameter-sweep gradient (r/d/b).
COLOR_TRUTH = "0.35"
COLOR_FLOW = "tab:blue"
COLOR_TRANSPORT = "tab:orange"
CMAP_SWEEP = "viridis"

PARAM_LABELS = {"discspacing": "discontinuity spacing $d$",
                "length_ratio": "length ratio $r$",
                "b": "aperture $b$"}

# Axis labels for the parity panels. The per-panel titles these replace were
# dropped: at readable type they collided, and they only repeated the labels.
MATH_LABEL = {"p0": r"$p_0$", "w1": r"$w$", "mu1": r"$\mu_1$",
              "sigma1": r"$\sigma_1$", "mu2": r"$\mu_2$", "sigma2": r"$\sigma_2$"}

SWEEP_ORDER = ["length_ratio", "discspacing", "b"]
N_SWEEP = 5


def load_transport_pipeline():
    """Re-run the transport notebook's data-loading + model-training
    prefix (cells 0-9, 71aaab88) in a fresh namespace, returning it as a
    dict of the variables the figures need."""
    nb_path = UPDATED_CODE_DIR / "[per p, t] probability of >=p percent having bkt <= t.ipynb"
    nb = json.load(open(nb_path))
    cells = {c["id"]: "".join(c["source"]) for c in nb["cells"]}
    order = ["cell-0", "cell-1", "cell-2", "cell-3", "cell-4", "cell-5",
             "cell-6", "cell-7", "cell-9", "71aaab88"]
    src = "\n\n".join(cells[cid] for cid in order)
    ns = {}
    exec(compile(src, str(nb_path), "exec"), ns)
    return ns


def load_flow_pipeline():
    """Load the flow model's saved results table + trained NN, recomputing
    the (deterministic, data-derived) normalization stats needed to feed
    new synthetic theta into the NN."""
    import tensorflow as tf
    from tensorflow import keras

    outputs_dir = Path("/Users/angelaliu/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Fracture Flow UROP/outputs")
    flow_results_df = pd.read_csv(outputs_dir / "flow_results_df.csv")
    # compile=False -- flow_model.keras was saved with a custom weighted_mse
    # loss function that isn't registered for deserialization, and we only
    # need forward inference (predict/__call__) here, not further training.
    flow_model = keras.models.load_model(outputs_dir / "flow_model.keras", compile=False)

    thetas = flow_results_df[["discspacing", "length_ratio", "b"]].to_numpy(dtype=float)
    outputs = flow_results_df[["p0_true", "w1_true", "mu1_true", "sigma1_true",
                                "mu2_true", "sigma2_true"]].to_numpy(dtype=float)
    # logwd (the 7th NN output) isn't in flow_results_df -- reconstruct its
    # mean/std is unnecessary since we only need thetas_means/stds and
    # outputs_means/stds for the first 6 output columns to denormalize
    # flow_model's predictions; the 7th column (logwd) is never consumed
    # downstream of the NN's raw prediction here, but the NN outputs 7
    # values, so pad with the true logwd distribution's own stats to keep
    # denormalization correct for the columns we DO use.
    logwd_col = np.log(flow_results_df.get("wd", pd.Series(np.nan, index=flow_results_df.index)))
    outputs_full = np.column_stack([outputs, logwd_col.to_numpy()])

    thetas_means, thetas_stds = thetas.mean(axis=0), thetas.std(axis=0)
    outputs_means, outputs_stds = np.nanmean(outputs_full, axis=0), np.nanstd(outputs_full, axis=0)

    return dict(flow_results_df=flow_results_df, flow_model=flow_model,
                thetas_means=thetas_means, thetas_stds=thetas_stds,
                outputs_means=outputs_means, outputs_stds=outputs_stds)


def flow_predict_mixture(flow_ns, disc, lr, b):
    """theta -> (p0, w1, mu1, sigma1, mu2, sigma2) in log(Q)-space, denormalized."""
    theta_norm = (np.array([disc, lr, b]) - flow_ns["thetas_means"]) / flow_ns["thetas_stds"]
    raw = flow_ns["flow_model"](theta_norm.reshape(1, -1)).numpy().flatten()
    raw = raw * flow_ns["outputs_stds"] + flow_ns["outputs_means"]
    p0, w1, mu1, sigma1, mu2, sigma2, _logwd = raw
    return p0, w1, mu1, abs(sigma1), mu2, abs(sigma2)


def make_fig1(tp, fp):
    """Held-out error distribution for each surrogate, side by side."""
    K = scale_fonts_for(9.0)
    flow_val = fp["flow_results_df"].loc[fp["flow_results_df"]["is_val"], "wd"]
    error_score_df = tp["error_score_df"]
    val_ps = tp["val_ps"]
    transport_val = error_score_df.loc[error_score_df["paramset_idx"].isin(val_ps),
                                       "error_score_constrained"]

    fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6))
    fig.set_layout_engine("constrained")

    axes[0].hist(np.log10(flow_val), bins=15, color=COLOR_FLOW, edgecolor="white")
    axes[0].axvline(np.log10(flow_val.median()), color="black", ls="--",
                    label=f"median = {flow_val.median():.2e}")
    axes[0].set_xlabel(r"$\log_{10}(\mathrm{WD})$")
    axes[0].set_ylabel("count")
    axes[0].set_title(f"Flow model (n={len(flow_val)})")
    axes[0].legend()

    axes[1].hist(np.log10(transport_val), bins=15, color=COLOR_TRANSPORT, edgecolor="white")
    axes[1].axvline(np.log10(transport_val.median()), color="black", ls="--",
                    label=f"median = {transport_val.median():.3f}")
    axes[1].set_xlabel(r"$\log_{10}(e)$")
    axes[1].set_ylabel("count")
    axes[1].set_title(f"Transport model (n={len(transport_val)})")
    axes[1].legend()

    fig.savefig(FIGURES_DIR / "fig1_aggregate_accuracy.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_aggregate_accuracy.png")


def make_fig2a_parity_flow(fp):
    """Predicted vs simulated mixture coefficients, held-out parameter sets."""
    K = scale_fonts_for(9.0)
    flow_val_df = fp["flow_results_df"][fp["flow_results_df"]["is_val"]]
    panels = ["p0", "w1", "mu1", "sigma1", "mu2", "sigma2"]

    fig, axes = plt.subplots(2, 3, figsize=(9.0, 5.6))
    fig.set_layout_engine("constrained")
    for name, ax in zip(panels, axes.ravel()):
        true_v, pred_v = flow_val_df[f"{name}_true"], flow_val_df[f"{name}_pred"]
        ax.scatter(true_v, pred_v, s=12 * K**2, alpha=0.7, color=COLOR_FLOW)
        lims = [min(true_v.min(), pred_v.min()), max(true_v.max(), pred_v.max())]
        ax.plot(lims, lims, "k--", lw=1 * K, alpha=0.6)
        ax.set_xlabel(f"true {MATH_LABEL[name]}")
        ax.set_ylabel(f"predicted {MATH_LABEL[name]}")

    fig.savefig(FIGURES_DIR / "fig2a_parity_flow.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("Saved fig2a_parity_flow.png")


def _transport_parity_points(tp):
    """Pooled (empirical F, predicted F) pairs over the held-out parameter sets."""
    true_vals, pred_vals = [], []
    for ps_idx in tp["val_ps"]:
        disc, lr, b = tp["recover_theta"](tp["weights"][ps_idx])
        for q in tp["q_vals"]:
            result = tp["group_ecdf"](ps_idx, q)
            if result is None:
                continue
            log_vals, ecdf = result
            X_g = tp["normalize_X"](np.column_stack([
                np.full(len(log_vals), disc), np.full(len(log_vals), lr),
                np.full(len(log_vals), b),
                np.full(len(log_vals), tp["logit"](q)), log_vals,
            ])).astype(np.float32)
            pred_vals.append(tp["model_constrained"](X_g, training=False).numpy().squeeze(axis=1))
            true_vals.append(ecdf)
    return np.concatenate(true_vals), np.concatenate(pred_vals)


def make_fig2b_parity_transport(tp):
    """Pooled predicted-vs-empirical F for the constrained transport model."""
    K = scale_fonts_for(5.0)
    true_vals, pred_vals = _transport_parity_points(tp)

    fig, ax = plt.subplots(figsize=(5.0, 4.6))
    fig.set_layout_engine("constrained")
    ax.scatter(true_vals, pred_vals, s=3 * K**2, alpha=0.15,
               color=COLOR_TRANSPORT, rasterized=True)
    ax.plot([0, 1], [0, 1], "k--", lw=1.2 * K, alpha=0.7)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_aspect("equal")
    ax.set_xlabel("empirical $F$")
    ax.set_ylabel("predicted $F$")

    fig.savefig(FIGURES_DIR / "fig2b_parity_transport.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print(f"Saved fig2b_parity_transport.png  ({len(true_vals)} points, "
          f"{len(tp['val_ps'])} held-out param sets)")


def make_fig6(tp, fp):
    """Held-out error against each theta component. Per-panel titles are
    omitted -- they duplicated the axis labels and collided once the type was
    large enough to read; the row identity now sits in the column-0 y-label."""
    K = scale_fonts_for(10.0)
    error_score_df = tp["error_score_df"]
    transport_val_df = error_score_df[error_score_df["paramset_idx"].isin(tp["val_ps"])]
    flow_val_df = fp["flow_results_df"][fp["flow_results_df"]["is_val"]]

    fig, axes = plt.subplots(2, 3, figsize=(10.0, 5.6), sharey="row")
    fig.set_layout_engine("constrained")
    for col, tcol in enumerate(["discspacing", "length_ratio", "b"]):
        ax = axes[0, col]
        ax.scatter(flow_val_df[tcol], flow_val_df["wd"], s=12 * K**2,
                   alpha=0.7, color=COLOR_FLOW)
        ax.set_yscale("log")
        ax.set_xlabel(PARAM_LABELS[tcol])
        if col == 0:
            ax.set_ylabel("flow model:\n" + r"$\mathrm{WD}$")

        ax2 = axes[1, col]
        ax2.scatter(transport_val_df[tcol], transport_val_df["error_score_constrained"],
                    s=12 * K**2, alpha=0.7, color=COLOR_TRANSPORT)
        ax2.set_yscale("log")
        ax2.set_xlabel(PARAM_LABELS[tcol])
        if col == 0:
            ax2.set_ylabel("transport model:\n" + r"error score $e$")

    fig.savefig(FIGURES_DIR / "fig6_error_vs_parameters.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig6_error_vs_parameters.png")


def _sweep_values(df, sweep_param, other_params):
    """Fixed values (training median) and sweep range (10th-90th percentile),
    taken from whichever model's own training distribution df describes."""
    fixed = {p: df[p].median() for p in other_params}
    sweep = np.linspace(df[sweep_param].quantile(0.1),
                        df[sweep_param].quantile(0.9), N_SWEEP)
    return fixed, sweep


def make_sweeps_flow(fp):
    """Flow-surrogate response to each parameter, as one 1x3 figure.

    Replaces the left-hand panels of the old per-parameter fig3/fig4/fig5."""
    K = scale_fonts_for(12.0)
    colors = matplotlib.colormaps[CMAP_SWEEP](np.linspace(0.05, 0.95, N_SWEEP))
    flow_df = fp["flow_results_df"]
    x_logq = np.linspace(np.log(1e-8), np.log(1e-2), 400)

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    fig.set_layout_engine("constrained")
    for ax, sweep_param in zip(axes, SWEEP_ORDER):
        others = [p for p in SWEEP_ORDER if p != sweep_param]
        fixed, sweep = _sweep_values(flow_df, sweep_param, others)
        for val, color in zip(sweep, colors):
            theta = dict(fixed, **{sweep_param: val})
            p0, w1, mu1, sigma1, mu2, sigma2 = flow_predict_mixture(
                fp, theta["discspacing"], theta["length_ratio"], theta["b"])
            pdf = (w1 * scipy_norm.pdf(x_logq, mu1, sigma1)
                   + (1 - w1) * scipy_norm.pdf(x_logq, mu2, sigma2))
            ax.plot(x_logq, pdf, color=color, lw=1.8 * K)
        ax.set_xlabel(r"$\log Q$")
        # a compact colorbar reads better here than a 5-entry legend
        sm = plt.cm.ScalarMappable(cmap=CMAP_SWEEP,
                                   norm=plt.Normalize(vmin=sweep.min(), vmax=sweep.max()))
        cbar = fig.colorbar(sm, ax=ax, orientation="horizontal",
                            location="top", pad=0.02, fraction=0.09)
        cbar.set_label(PARAM_LABELS[sweep_param])
        cbar.ax.tick_params(labelsize=TARGET_PT["xtick.labelsize"] * K)
    axes[0].set_ylabel("predicted PDF")

    fig.savefig(FIGURES_DIR / "fig7_sweeps_flow.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("Saved fig7_sweeps_flow.png")


def make_sweeps_transport(tp):
    """Transport-surrogate response to each parameter, as one 1x3 figure.

    Replaces the right-hand panels of the old per-parameter fig3/fig4/fig5."""
    K = scale_fonts_for(12.0)
    colors = matplotlib.colormaps[CMAP_SWEEP](np.linspace(0.05, 0.95, N_SWEEP))

    thetas_usable = np.array([tp["recover_theta"](tp["weights"][i]) for i in tp["USABLE_PS"]])
    df_theta = pd.DataFrame(thetas_usable, columns=["discspacing", "length_ratio", "b"])
    T_GRID = np.linspace(tp["flat_df"]["logt"].min(), tp["flat_df"]["logt"].max(), 60)
    Q_GRID = tp["Q_GRID"]
    Tm, Qm = np.meshgrid(T_GRID, Q_GRID)
    predict_fn = lambda X: tp["model_constrained"](X, training=False).numpy().squeeze(axis=1)

    fig = plt.figure(figsize=(12.0, 4.4))
    for i, sweep_param in enumerate(SWEEP_ORDER):
        ax = fig.add_subplot(1, 3, i + 1, projection="3d")
        others = [p for p in SWEEP_ORDER if p != sweep_param]
        fixed, sweep = _sweep_values(df_theta, sweep_param, others)
        for val, color in zip(sweep, colors):
            theta = dict(fixed, **{sweep_param: val})
            X_grid = tp["normalize_X"](np.column_stack([
                np.full(Tm.size, theta["discspacing"]),
                np.full(Tm.size, theta["length_ratio"]),
                np.full(Tm.size, theta["b"]),
                tp["logit"](Qm.ravel()), Tm.ravel(),
            ])).astype(np.float32)
            Z = predict_fn(X_grid).reshape(Tm.shape)
            ax.plot_surface(Tm, Qm, Z, color=color, alpha=0.55, linewidth=0, antialiased=True)
        ax.set_xlabel(r"$\log t$", labelpad=6 * K)
        ax.set_ylabel(r"$q$", labelpad=4 * K)
        ax.set_zlabel(r"$F$", labelpad=3 * K)
        ax.set_zlim(0, 1)
        ax.view_init(elev=22, azim=-55)
        ax.tick_params(pad=1.0 * K)
        # fewer ticks: at this panel size the log-t and q tick labels collide
        # at the shared corner of the 3D box
        ax.xaxis.set_major_locator(matplotlib.ticker.MaxNLocator(4))
        ax.set_yticks([0.0, 0.5, 1.0])
        ax.set_zticks([0.0, 0.5, 1.0])
        # shrink the 3D content inside its axes so the z-label is not clipped
        # by bbox_inches="tight", which crops to the axes box
        ax.set_box_aspect(None, zoom=0.86)

        sm = plt.cm.ScalarMappable(cmap=CMAP_SWEEP,
                                   norm=plt.Normalize(vmin=sweep.min(), vmax=sweep.max()))
        cbar = fig.colorbar(sm, ax=ax, orientation="horizontal",
                            location="top", pad=0.0, fraction=0.05, shrink=0.85)
        cbar.set_label(PARAM_LABELS[sweep_param])
        cbar.ax.tick_params(labelsize=TARGET_PT["xtick.labelsize"] * K)

    fig.subplots_adjust(left=0.01, right=0.99, bottom=0.02, top=0.90, wspace=0.22)
    fig.savefig(FIGURES_DIR / "fig8_sweeps_transport.png", bbox_inches="tight", dpi=200)
    plt.close(fig)
    print("Saved fig8_sweeps_transport.png")


def make_lambda_sweep():
    """Heatmaps of the transport model's 2D physics-constraint-weight sweep
    (lam_t, lam_q) -- reads the cached CSV written by the notebook's
    lamsweep01 cell (25 retrains, ~8-10 min) rather than re-running the
    sweep here, so this figure regenerates in seconds.

    Renamed from fig7_lambda_sweep.png to figA1_lambda_sweep.png: fig7 now
    belongs to the flow sweeps, and this figure is no longer cited in
    Section 4 (it supports the appendix / the lambda choice in Section 3)."""
    K = scale_fonts_for(11.0)
    transport_dir = Path("/Users/angelaliu/Library/CloudStorage/OneDrive-MassachusettsInstituteofTechnology/Fracture Flow UROP/Updated Code/transport")
    csv_path = transport_dir / "lam_sweep_results.csv"
    if not csv_path.exists():
        print(f"Skipping figA1_lambda_sweep.png -- {csv_path} not found "
              f"(run the notebook's lamsweep01 cell at least once to produce it).")
        return

    sweep_df = pd.read_csv(csv_path)
    e_grid = sweep_df.pivot(index="lam_t", columns="lam_q", values="e").values
    viol_grid = sweep_df.pivot(index="lam_t", columns="lam_q", values="q_viol_rate").values * 100
    lam_t_vals = sorted(sweep_df["lam_t"].unique())
    lam_q_vals = sorted(sweep_df["lam_q"].unique())
    best_row = sweep_df.loc[sweep_df["e"].idxmin()]

    fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.4))
    fig.set_layout_engine("constrained")
    for ax, grid, title, cbar_label, fmt in zip(
        axes, [e_grid, viol_grid],
        [r"mean $e$", r"$\partial F/\partial q>0$ violation rate"],
        ["mean $e$", "violation rate (%)"],
        ["{:.3f}", "{:.2f}%"],
    ):
        im = ax.imshow(grid, cmap=CMAP_SWEEP, aspect="auto")
        ax.set_xticks(range(len(lam_q_vals)), lam_q_vals)
        ax.set_yticks(range(len(lam_t_vals)), lam_t_vals)
        ax.set_xlabel(r"$\lambda_q$")
        ax.set_ylabel(r"$\lambda_t$")
        ax.set_title(title)
        thr = grid.min() + 0.5 * (grid.max() - grid.min())
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                ax.text(j, i, fmt.format(grid[i, j]), ha="center", va="center",
                        color="black" if grid[i, j] > thr else "white",
                        fontsize=TARGET_PT["xtick.labelsize"] * K)
        cbar = fig.colorbar(im, ax=ax)
        cbar.set_label(cbar_label)
        cbar.ax.tick_params(labelsize=TARGET_PT["ytick.labelsize"] * K)

    fig.suptitle(f"Constraint-weight sweep (best: $\\lambda_t$={best_row['lam_t']:g}, "
                 f"$\\lambda_q$={best_row['lam_q']:g}, $e$={best_row['e']:.3f})")
    fig.savefig(FIGURES_DIR / "figA1_lambda_sweep.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved figA1_lambda_sweep.png")


if __name__ == "__main__":
    print("Loading transport-model pipeline (retraining model_plain/model_constrained)...")
    tp = load_transport_pipeline()
    print("Loading flow-model pipeline (saved results + NN)...")
    fp = load_flow_pipeline()
    print("Both pipelines loaded. Generating figures...")

    make_fig1(tp, fp)
    make_fig2a_parity_flow(fp)
    make_fig2b_parity_transport(tp)
    make_fig6(tp, fp)
    make_sweeps_flow(fp)
    make_sweeps_transport(tp)
    make_lambda_sweep()
    print("All figures generated.")
