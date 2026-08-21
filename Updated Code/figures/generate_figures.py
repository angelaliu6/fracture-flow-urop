"""Generate the summary figure set (fig1-fig6) comparing the flow-model
surrogate ("Updated Code/flow distrib.ipynb") and the transport-model
surrogate ("Updated Code/[per p, t] probability of >=p percent having bkt
<= t.ipynb"). Run standalone: `.venv/bin/python "Updated Code/figures/generate_figures.py"`.

The transport-model pipeline (model_constrained, error_score_df, val_ps,
and the plotting helpers) isn't saved as a loadable artifact anywhere --
only the .keras weights are -- so it's regenerated here by re-running the
relevant prefix of that notebook's own cells in-process (same values each
time: no randomness in the data/normalization, and retraining the NNs
takes ~15-20s and gives an equivalent-quality model on every run). This
keeps a single source of truth (the notebook) instead of a second,
driftable copy of that logic living in this script.
"""
import json
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from pathlib import Path

matplotlib.rcParams.update({
    "font.size": 14,
    "axes.titlesize": 17,
    "axes.labelsize": 15,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.dpi": 150,
})

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
    flow_val = fp["flow_results_df"].loc[fp["flow_results_df"]["is_val"], "wd"]
    error_score_df = tp["error_score_df"]
    val_ps = tp["val_ps"]
    transport_val = error_score_df.loc[error_score_df["paramset_idx"].isin(val_ps), "error_score_constrained"]

    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    axes[0].hist(np.log10(flow_val), bins=15, color=COLOR_FLOW, edgecolor="white")
    axes[0].axvline(np.log10(flow_val.median()), color="black", ls="--", lw=1.5,
                     label=f"median = {flow_val.median():.2e}")
    axes[0].set_xlabel(r"$\log_{10}(w_d)$  (flow model, held-out)")
    axes[0].set_ylabel("count")
    axes[0].set_title(f"Flow model: held-out Wasserstein distance\n(n={len(flow_val)} param sets)")
    axes[0].legend()

    axes[1].hist(np.log10(transport_val), bins=15, color=COLOR_TRANSPORT, edgecolor="white")
    axes[1].axvline(np.log10(transport_val.median()), color="black", ls="--", lw=1.5,
                     label=f"median = {transport_val.median():.4f}")
    axes[1].set_xlabel(r"$\log_{10}(e)$  (transport model, held-out)")
    axes[1].set_ylabel("count")
    axes[1].set_title(f"Transport model: held-out error score\n(n={len(transport_val)} param sets)")
    axes[1].legend()

    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig1_aggregate_accuracy.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig1_aggregate_accuracy.png")


def make_fig2(tp, fp):
    fig = plt.figure(figsize=(20, 9))
    gs = fig.add_gridspec(2, 5, width_ratios=[1, 1, 1, 0.15, 2])

    flow_val_df = fp["flow_results_df"][fp["flow_results_df"]["is_val"]]
    flow_panels = [("p0", "p0_true", "p0_pred"), ("w1", "w1_true", "w1_pred"),
                    ("mu1", "mu1_true", "mu1_pred"), ("sigma1", "sigma1_true", "sigma1_pred"),
                    ("mu2", "mu2_true", "mu2_pred"), ("sigma2", "sigma2_true", "sigma2_pred")]
    positions = [(0, 0), (0, 1), (0, 2), (1, 0), (1, 1), (1, 2)]
    for (name, true_col, pred_col), (r, c) in zip(flow_panels, positions):
        ax = fig.add_subplot(gs[r, c])
        true_v, pred_v = flow_val_df[true_col], flow_val_df[pred_col]
        ax.scatter(true_v, pred_v, s=22, alpha=0.7, color=COLOR_FLOW)
        lims = [min(true_v.min(), pred_v.min()), max(true_v.max(), pred_v.max())]
        ax.plot(lims, lims, "k--", lw=1, alpha=0.6)
        ax.set_xlabel(f"true {name}")
        ax.set_ylabel(f"predicted {name}")
        ax.set_title(name, fontsize=13)

    ax_transport = fig.add_subplot(gs[:, 4])
    error_pt_df = tp["error_pt_df"]
    val_ps = tp["val_ps"]
    true_vals, pred_vals = [], []
    for ps_idx in val_ps:
        disc, lr, b = tp["recover_theta"](tp["weights"][ps_idx])
        for p in tp["p_vals"]:
            result = tp["group_ecdf"](ps_idx, p)
            if result is None:
                continue
            log_vals, ecdf = result
            X_g = tp["normalize_X"](np.column_stack([
                np.full(len(log_vals), disc), np.full(len(log_vals), lr), np.full(len(log_vals), b),
                np.full(len(log_vals), tp["logit"](p)), log_vals,
            ])).astype(np.float32)
            F_pred = tp["model_constrained"](X_g, training=False).numpy().squeeze(axis=1)
            true_vals.append(ecdf)
            pred_vals.append(F_pred)
    true_vals = np.concatenate(true_vals)
    pred_vals = np.concatenate(pred_vals)
    ax_transport.scatter(true_vals, pred_vals, s=3, alpha=0.15, color=COLOR_TRANSPORT, rasterized=True)
    ax_transport.plot([0, 1], [0, 1], "k--", lw=1.5, alpha=0.7)
    ax_transport.set_xlim(0, 1)
    ax_transport.set_ylim(0, 1)
    ax_transport.set_xlabel("actual $F$ (empirical ecdf)")
    ax_transport.set_ylabel("predicted $F$ (model_constrained)")
    ax_transport.set_title(f"Transport model: pooled parity\n({len(true_vals)} (theta,q,t) points, {len(val_ps)} held-out param sets)")

    fig.suptitle("Parity: predicted vs. true, held-out param sets", fontsize=19)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig2_parity.png", bbox_inches="tight", dpi=130)
    plt.close(fig)
    print("Saved fig2_parity.png")


from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from scipy.stats import norm as scipy_norm


PARAM_LABELS = {"discspacing": "discspacing $d$", "length_ratio": "length ratio $r$", "b": "aperture $b$"}
PARAM_FIG_NAMES = {"discspacing": "fig4_impact_of_discspacing.png",
                    "length_ratio": "fig3_impact_of_length_ratio.png",
                    "b": "fig5_impact_of_b.png"}


def impact_of_parameter(tp, fp, sweep_param):
    """Two-panel figure at fixed representative values of the OTHER two
    theta components: (a) family of flow-model log(Q) PDF curves swept
    over sweep_param, viridis gradient; (b) family of transport-model 3D
    F(q,logt) surfaces over the same sweep, viridis gradient. Each panel
    uses ITS OWN dataset's representative fixed values and sweep range
    (10th-90th percentile of that param, from that model's own data) --
    the two datasets are independent samples over a similar but not
    identical range, so this keeps each sweep within its own model's
    training distribution."""
    other_params = [p for p in ["discspacing", "length_ratio", "b"] if p != sweep_param]
    N_SWEEP = 5
    sweep_colors = matplotlib.colormaps[CMAP_SWEEP](np.linspace(0.05, 0.95, N_SWEEP))

    fig = plt.figure(figsize=(18, 8))
    ax_flow = fig.add_subplot(1, 2, 1)
    ax_transport = fig.add_subplot(1, 2, 2, projection="3d")

    # --- (a) flow model ---
    flow_df = fp["flow_results_df"]
    fixed_flow = {p: flow_df[p].median() for p in other_params}
    sweep_flow = np.linspace(flow_df[sweep_param].quantile(0.1), flow_df[sweep_param].quantile(0.9), N_SWEEP)
    x_logq = np.linspace(np.log(1e-8), np.log(1e-2), 400)
    for val, color in zip(sweep_flow, sweep_colors):
        theta = dict(fixed_flow, **{sweep_param: val})
        p0, w1, mu1, sigma1, mu2, sigma2 = flow_predict_mixture(fp, theta["discspacing"], theta["length_ratio"], theta["b"])
        pdf = w1 * scipy_norm.pdf(x_logq, mu1, sigma1) + (1 - w1) * scipy_norm.pdf(x_logq, mu2, sigma2)
        ax_flow.plot(x_logq, pdf, color=color, lw=2, label=f"{val:.3g}")
    ax_flow.set_xlabel(r"$\log(Q)$  (network flow rate)")
    ax_flow.set_ylabel("predicted PDF")
    ax_flow.set_title(f"Flow model: predicted log(Q) distribution\n(fixed {', '.join(f'{k}={v:.3g}' for k,v in fixed_flow.items())})")
    ax_flow.legend(title=PARAM_LABELS[sweep_param], fontsize=10)

    # --- (b) transport model ---
    val_ps = tp["val_ps"]
    thetas_usable = np.array([tp["recover_theta"](tp["weights"][i]) for i in tp["USABLE_PS"]])
    df_theta = pd.DataFrame(thetas_usable, columns=["discspacing", "length_ratio", "b"])
    fixed_transport = {p: df_theta[p].median() for p in other_params}
    sweep_transport = np.linspace(df_theta[sweep_param].quantile(0.1), df_theta[sweep_param].quantile(0.9), N_SWEEP)
    T_GRID = np.linspace(tp["flat_df"]["logt"].min(), tp["flat_df"]["logt"].max(), 60)
    Q_GRID = tp["P_GRID"]
    Tm, Qm = np.meshgrid(T_GRID, Q_GRID)
    predict_fn = lambda X: tp["model_constrained"](X, training=False).numpy().squeeze(axis=1)
    for val, color in zip(sweep_transport, sweep_colors):
        theta = dict(fixed_transport, **{sweep_param: val})
        X_grid = tp["normalize_X"](np.column_stack([
            np.full(Tm.size, theta["discspacing"]), np.full(Tm.size, theta["length_ratio"]), np.full(Tm.size, theta["b"]),
            tp["logit"](Qm.ravel()), Tm.ravel(),
        ])).astype(np.float32)
        Z = predict_fn(X_grid).reshape(Tm.shape)
        ax_transport.plot_surface(Tm, Qm, Z, color=color, alpha=0.55, linewidth=0, antialiased=True)
    ax_transport.set_xlabel("log t", labelpad=10)
    ax_transport.set_ylabel("percentile $q$", labelpad=10)
    ax_transport.set_zlabel("predicted $F$", labelpad=6)
    ax_transport.set_zlim(0, 1)
    ax_transport.view_init(elev=22, azim=-55)
    ax_transport.set_title(f"Transport model: predicted $F(q,\\log t)$ surface\n(fixed {', '.join(f'{k}={v:.3g}' for k,v in fixed_transport.items())})")

    sm = plt.cm.ScalarMappable(cmap=CMAP_SWEEP, norm=plt.Normalize(vmin=sweep_transport.min(), vmax=sweep_transport.max()))
    cbar = fig.colorbar(sm, ax=ax_transport, shrink=0.6, pad=0.1)
    cbar.set_label(PARAM_LABELS[sweep_param])

    fig.suptitle(f"Impact of {PARAM_LABELS[sweep_param]}", fontsize=19)
    plt.tight_layout()
    fig_path = FIGURES_DIR / PARAM_FIG_NAMES[sweep_param]
    fig.savefig(fig_path, bbox_inches="tight", dpi=130)
    plt.close(fig)
    print(f"Saved {fig_path.name}")


def make_fig6(tp, fp):
    error_score_df = tp["error_score_df"]
    val_ps = tp["val_ps"]
    transport_val_df = error_score_df[error_score_df["paramset_idx"].isin(val_ps)]
    flow_val_df = fp["flow_results_df"][fp["flow_results_df"]["is_val"]]

    theta_cols = ["discspacing", "length_ratio", "b"]
    fig, axes = plt.subplots(2, 3, figsize=(17, 10))
    for col, tcol in enumerate(theta_cols):
        ax = axes[0, col]
        ax.scatter(flow_val_df[tcol], flow_val_df["wd"], s=22, alpha=0.7, color=COLOR_FLOW)
        ax.set_yscale("log")
        ax.set_xlabel(PARAM_LABELS[tcol])
        ax.set_ylabel(r"$w_d$" if col == 0 else "")
        ax.set_title(f"flow model: $w_d$ vs {PARAM_LABELS[tcol]}", fontsize=12)

        ax2 = axes[1, col]
        ax2.scatter(transport_val_df[tcol], transport_val_df["error_score_constrained"], s=22, alpha=0.7, color=COLOR_TRANSPORT)
        ax2.set_yscale("log")
        ax2.set_xlabel(PARAM_LABELS[tcol])
        ax2.set_ylabel("error score $e$" if col == 0 else "")
        ax2.set_title(f"transport model: $e$ vs {PARAM_LABELS[tcol]}", fontsize=12)

    fig.suptitle("Held-out prediction error vs. each theta component", fontsize=19)
    plt.tight_layout()
    fig.savefig(FIGURES_DIR / "fig6_error_vs_parameters.png", bbox_inches="tight")
    plt.close(fig)
    print("Saved fig6_error_vs_parameters.png")


if __name__ == "__main__":
    print("Loading transport-model pipeline (retraining model_plain/model_constrained)...")
    tp = load_transport_pipeline()
    print("Loading flow-model pipeline (saved results + NN)...")
    fp = load_flow_pipeline()
    print("Both pipelines loaded. Generating figures...")

    make_fig1(tp, fp)
    make_fig2(tp, fp)
    impact_of_parameter(tp, fp, "length_ratio")
    impact_of_parameter(tp, fp, "discspacing")
    impact_of_parameter(tp, fp, "b")
    make_fig6(tp, fp)
    print("All figures generated.")
