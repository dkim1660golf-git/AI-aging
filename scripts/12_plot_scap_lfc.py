from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

scap_deg_path = Path("results/scap_starter_deg_GSE346632.csv")
out_png = Path("results/scap_starter_log2fc.png")

df = pd.read_csv(scap_deg_path)

# log2FoldChange 기준 내림차순 정렬 (가장 큰 양수가 맨 위)
df = df.sort_values("log2FoldChange", ascending=False).reset_index(drop=True)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False
sns.set_theme(style="whitegrid", context="notebook")

colors = []
for _, row in df.iterrows():
    is_sig = pd.notna(row["padj"]) and (row["padj"] < 0.05)
    if not is_sig:
        colors.append("#94A3B8")  # NS
    elif row["log2FoldChange"] > 0:
        colors.append("#D95F02")  # UP (padj < 0.05, LFC > 0)
    else:
        colors.append("#2B6CB0")  # DOWN (padj < 0.05, LFC < 0)

fig, ax = plt.subplots(figsize=(9.5, 5.8))
fig.patch.set_facecolor("white")

y_pos = np.arange(len(df))
lfc_vals = df["log2FoldChange"].values
labels = [
    f"{r['symbol']} ({r['pathway_group']})" if "pathway_group" in df.columns and pd.notna(r.get("pathway_group")) else r["symbol"]
    for _, r in df.iterrows()
]

bars = ax.barh(y_pos, lfc_vals, height=0.62, color=colors, edgecolor="none")
ax.axvline(0.0, color="#1E293B", linewidth=1.2, linestyle="-")

# 각 막대 끝에 log2FC 값과 유의성 표시
x_min, x_max = lfc_vals.min(), lfc_vals.max()
span = max(abs(x_min), abs(x_max))
pad = span * 0.03

for i, (_, row) in enumerate(df.iterrows()):
    val = float(row["log2FoldChange"])
    is_sig = pd.notna(row["padj"]) and (row["padj"] < 0.05)
    star = "*" if is_sig else " (ns)"
    txt = f"{val:+.2f}{star}"

    if val >= 0:
        ax.text(
            val + pad,
            i,
            txt,
            va="center",
            ha="left",
            fontsize=9,
            fontweight="bold" if is_sig else "normal",
            color="#1E293B" if is_sig else "#64748B",
        )
    else:
        ax.text(
            val - pad,
            i,
            txt,
            va="center",
            ha="right",
            fontsize=9,
            fontweight="bold" if is_sig else "normal",
            color="#1E293B" if is_sig else "#64748B",
        )

ax.set_yticks(y_pos)
ax.set_yticklabels(labels, fontsize=9.5, fontweight="bold")
ax.invert_yaxis()

ax.set_xlim(x_min - span * 0.25, x_max + span * 0.25)
ax.set_xlabel("log2FoldChange (PD40 vs PD32)", fontsize=10.5, fontweight="bold")
ax.set_title(
    "SCAP Starter Genes in GSE346632 (IMR90: PD40 vs PD32)",
    fontsize=12.5,
    fontweight="bold",
    loc="left",
    pad=12,
)

# 범례 추가
legend_handles = [
    plt.Rectangle((0, 0), 1, 1, facecolor="#D95F02", label="UP (padj < 0.05, LFC > 0)"),
    plt.Rectangle((0, 0), 1, 1, facecolor="#2B6CB0", label="DOWN (padj < 0.05, LFC < 0)"),
    plt.Rectangle((0, 0), 1, 1, facecolor="#94A3B8", label="Not Significant (ns)"),
]
ax.legend(
    handles=legend_handles,
    loc="lower right",
    frameon=True,
    facecolor="white",
    edgecolor="#E2E8F0",
    fontsize=8.5,
)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()

print("SCAP log2FC 그래프 저장 완료:", out_png)