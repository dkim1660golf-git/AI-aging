from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

summary_csv = Path("results/l1000fwd/02_drug_summary_chembl.csv")
fallback_csv = Path("results/l1000fwd/04_final_repurposing_candidates.csv")
out_png = Path("results/l1000fwd/07_real_reversal_significance_map.png")

src_csv = summary_csv if summary_csv.exists() else fallback_csv
if not src_csv.exists():
    raise FileNotFoundError("시각화할 CSV 파일을 찾지 못했습니다.")

df = pd.read_csv(src_csv)
if df.empty:
    raise ValueError("CSV 데이터가 비어 있습니다.")

df["min_similarity_score"] = pd.to_numeric(df.get("min_similarity_score", 0), errors="coerce").fillna(0.0)
df["best_q_value"] = pd.to_numeric(df.get("best_q_value", np.nan), errors="coerce")
df["max_phase_num"] = pd.to_numeric(df.get("max_phase", 0), errors="coerce").fillna(0.0)
df["abs_sim"] = df["min_similarity_score"].abs()

valid_q = df["best_q_value"].where(df["best_q_value"] > 0).clip(lower=1e-12, upper=1.0).fillna(1.0)
df["neg_log10_q"] = -np.log10(valid_q)

# 실제 similarity_score나 q-value 신호가 있는 행 우선 선택
signal_df = df[(df["abs_sim"] > 0) | (df["neg_log10_q"] > 0)].copy()
if signal_df.empty:
    signal_df = df.copy()

plot_df = (
    signal_df.sort_values(
        ["neg_log10_q", "abs_sim", "max_phase_num"],
        ascending=[False, False, False],
    )
    .head(12)
    .reset_index(drop=True)
)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

fig_h = max(4.8, 0.52 * len(plot_df) + 1.6)
fig, ax = plt.subplots(figsize=(11.5, fig_h))
fig.patch.set_facecolor("white")

y_pos = np.arange(len(plot_df))
x_vals = plot_df["neg_log10_q"].values
bubble_sizes = (plot_df["abs_sim"] * 12000.0).clip(lower=120.0, upper=900.0)
colors = ["#2B6CB0" if p >= 4 else ("#D95F02" if p > 0 else "#64748B") for p in plot_df["max_phase_num"]]

ax.scatter(
    x_vals,
    y_pos,
    s=bubble_sizes,
    c=colors,
    alpha=0.88,
    edgecolors="#1E293B",
    linewidth=1.2,
    zorder=3,
)

q_cutoff = -np.log10(0.05)
ax.axvline(q_cutoff, color="#DC2626", linestyle="--", linewidth=1.2, label="q = 0.05 significance cutoff", zorder=2)

x_max = max(float(np.nanmax(x_vals)) if len(x_vals) > 0 else 0.0, q_cutoff * 1.15)
x_right = x_max * 1.35 if np.isfinite(x_max) and x_max > 0 else 3.0
offset = x_right * 0.02

y_labels = []
for i, row in plot_df.iterrows():
    phase_str = f"Phase {int(row['max_phase_num'])}" if row["max_phase_num"] > 0 else "Unannotated"
    y_labels.append(f"{row['drug']} [{phase_str}]")

    moa_txt = str(row.get("moa", "")).strip()
    if not moa_txt or moa_txt.lower() == "nan":
        moa_txt = "No MOA recorded"
    else:
        moa_txt = moa_txt.split(";")[0].strip()
        if len(moa_txt) > 36:
            moa_txt = moa_txt[:33].rstrip() + "..."

    ax.text(
        row["neg_log10_q"] + offset,
        i,
        f"sim={row['min_similarity_score']:.4f} | {moa_txt}",
        va="center",
        fontsize=8.6,
        fontweight="bold",
        color="#1E293B",
    )

ax.set_yticks(y_pos)
ax.set_yticklabels(y_labels, fontsize=9.2, fontweight="bold")
ax.invert_yaxis()
ax.set_xlim(0, x_right)
ax.set_xlabel("-log10(best q-value)  [Bubble size = |Opposite Similarity Score|]", fontsize=9.8, fontweight="bold")
ax.set_title(
    "Real L1000FWD Opposite Signature Compounds (GSE346632: PD40 vs PD32)",
    fontsize=12,
    fontweight="bold",
    loc="left",
    pad=12,
)

legend_handles = [
    plt.Line2D([0], [0], marker="o", color="w", label="Approved (Phase 4)", markerfacecolor="#2B6CB0", markeredgecolor="#1E293B", markersize=9.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Clinical Trial (Phase 1-3)", markerfacecolor="#D95F02", markeredgecolor="#1E293B", markersize=9.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Preclinical / Unannotated", markerfacecolor="#64748B", markeredgecolor="#1E293B", markersize=9.5),
    plt.Line2D([0], [0], color="#DC2626", linestyle="--", linewidth=1.2, label="q = 0.05 cutoff"),
]
ax.legend(handles=legend_handles, loc="lower right", frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=8.3)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linestyle=":", alpha=0.5, zorder=0)

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()

print("실제 데이터 기준 반전 후보 맵 저장 완료:", out_png)