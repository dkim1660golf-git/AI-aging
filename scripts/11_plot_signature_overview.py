from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

summary_path = Path("results/signature_deg_summary.csv")
sample_scores_path = Path("results/sample_signature_scores.csv")
out_png = Path("results/signature_overview.png")

summary_df = pd.read_csv(summary_path)
scores_df = pd.read_csv(sample_scores_path)

# 깔끔한 전역 스타일 설정
plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
sns.set_theme(style="whitegrid", context="notebook")

sig_order = ["SenMayo", "SCAP_starter", "CellAge", "SASP_Atlas"]
summary_df["sig_rank"] = summary_df["signature"].apply(
    lambda x: sig_order.index(x) if x in sig_order else 99
)
summary_df = summary_df.sort_values("sig_rank").reset_index(drop=True)

sample_order = [
    "PD32-REP1",
    "PD32-REP2",
    "PD32-REP3",
    "PD40-REP1",
    "PD40-REP2",
    "PD40-REP3",
]

pivot_z = (
    scores_df.pivot(index="signature", columns="sample", values="mean_zscore")
    .reindex(index=sig_order, columns=sample_order)
)

fig, axes = plt.subplots(
    1, 2, figsize=(13.5, 5.2), gridspec_kw={"width_ratios": [1.05, 1.25]}
)
fig.patch.set_facecolor("white")

# [Panel A] 유의미한 상향(UP) vs 하향(DOWN) 유전자 수
y_pos = np.arange(len(summary_df))
up_vals = summary_df["up_sig_genes"].values
down_vals = summary_df["down_sig_genes"].values
labels = summary_df["signature"].values

axes[0].barh(
    y_pos - 0.17,
    up_vals,
    height=0.34,
    color="#D95F02",
    edgecolor="none",
    label="UP (padj<0.05, LFC>0)",
)
axes[0].barh(
    y_pos + 0.17,
    down_vals,
    height=0.34,
    color="#2B6CB0",
    edgecolor="none",
    label="DOWN (padj<0.05, LFC<0)",
)

for i, (u, d) in enumerate(zip(up_vals, down_vals)):
    axes[0].text(u + max(up_vals.max(), down_vals.max()) * 0.015, i - 0.17, f"{u}", va="center", fontsize=9, color="#9C4200", fontweight="bold")
    axes[0].text(d + max(up_vals.max(), down_vals.max()) * 0.015, i + 0.17, f"{d}", va="center", fontsize=9, color="#1E4E8C", fontweight="bold")

axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(labels, fontsize=10, fontweight="bold")
axes[0].invert_yaxis()
axes[0].set_xlabel("Significant DEG Count (padj < 0.05)", fontsize=10)
axes[0].set_title("A. Significant DEGs per Signature", fontsize=11.5, fontweight="bold", loc="left", pad=10)
axes[0].spines["top"].set_visible(False)
axes[0].spines["right"].set_visible(False)
axes[0].legend(frameon=True, facecolor="white", edgecolor="#E2E8F0", fontsize=8.5, loc="lower right")

# [Panel B] 샘플별 mean z-score 히트맵
hm = sns.heatmap(
    pivot_z,
    ax=axes[1],
    cmap="RdBu_r",
    center=0.0,
    annot=True,
    fmt=".2f",
    annot_kws={"size": 9, "weight": "bold"},
    linewidths=1.2,
    linecolor="white",
    cbar_kws={"label": "Mean Z-score", "shrink": 0.85},
)
axes[1].set_title("B. Sample Signature Z-score (PD32 vs PD40)", fontsize=11.5, fontweight="bold", loc="left", pad=10)
axes[1].set_xlabel("")
axes[1].set_ylabel("")
axes[1].tick_params(axis="x", rotation=25, labelsize=9.5)
axes[1].tick_params(axis="y", rotation=0, labelsize=10)

# PD32와 PD40 그룹 경계선 강조
axes[1].axvline(3.0, color="#1A202C", linestyle="--", linewidth=1.5)

fig.suptitle(
    "GSE346632 (IMR90 Replicative Senescence: PD40 vs PD32) Signature Overview",
    fontsize=13.5,
    fontweight="bold",
    y=0.99,
)
plt.tight_layout(rect=[0, 0, 1, 0.94])
plt.savefig(out_png, dpi=300)
plt.close()

print("그래프 저장 완료:", out_png)