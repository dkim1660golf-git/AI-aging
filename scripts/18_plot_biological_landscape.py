from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

in_csv = Path("results/l1000fwd/04_final_repurposing_candidates.csv")
out_png = Path("results/l1000fwd/06_biological_repurposing_landscape.png")

if not in_csv.exists():
    raise FileNotFoundError(f"입력 파일이 없습니다: {in_csv}")

df = pd.read_csv(in_csv)
if df.empty:
    raise ValueError("후보 데이터가 비어 있습니다.")

df["min_similarity_score"] = pd.to_numeric(df["min_similarity_score"], errors="coerce").fillna(0.0)
df["best_q_value"] = pd.to_numeric(df["best_q_value"], errors="coerce")
df["max_phase_num"] = pd.to_numeric(df.get("max_phase", 0), errors="coerce").fillna(0.0)
df["abs_sim"] = df["min_similarity_score"].abs()

valid_q = df["best_q_value"].where(df["best_q_value"] > 0).clip(lower=1e-12, upper=1.0)
df["neg_log10_q"] = -np.log10(valid_q)


def classify_bio_axis(moa_text: str) -> str:
    txt = str(moa_text).lower()
    if any(k in txt for k in ["tubulin", "microtubule", "kinesin", "cytoskeleton"]):
        return "Cytoskeleton / Mitosis"
    if any(k in txt for k in ["hsp", "heat shock", "proteasome", "chaperone"]):
        return "Proteostasis / HSP90"
    if any(k in txt for k in ["jak", "stat", "interleukin", "nf-kb", "cyclooxygenase", "cox"]):
        return "Inflammation / SASP Axis"
    if any(k in txt for k in ["pi3k", "akt", "mtor", "egfr", "erbb", "kinase", "raf", "mek"]):
        return "Survival / Kinase Signaling"
    if any(k in txt for k in ["retinoid", "retinoic", "vitamin d", "ppar", "histone", "hdac"]):
        return "Epigenetic / Metabolic Aging"
    return "Other / Unclassified MOA"


df["bio_axis"] = df["moa"].apply(classify_bio_axis) if "moa" in df.columns else "Other / Unclassified MOA"

# 점수가 있는 행을 우선 정렬하되, 모두 0/NaN이면 max_phase와 약물명 순으로 정렬
has_signal = (df["abs_sim"] > 0) | (df["neg_log10_q"].notna())
if has_signal.any():
    plot_df = df[has_signal].copy()
    plot_df = plot_df.sort_values(
        ["neg_log10_q", "abs_sim", "max_phase_num"],
        ascending=[False, False, False],
        na_position="last",
    ).head(10).reset_index(drop=True)
else:
    plot_df = df.sort_values(["max_phase_num", "drug"], ascending=[False, True]).head(10).reset_index(drop=True)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

fig_h = max(4.8, 0.55 * len(plot_df) + 1.6)
fig, ax = plt.subplots(figsize=(11.5, fig_h))
fig.patch.set_facecolor("white")

y_pos = np.arange(len(plot_df))
x_vals = plot_df["abs_sim"].values
point_colors = ["#2B6CB0" if p >= 4 else "#D95F02" for p in plot_df["max_phase_num"]]

q_sizes = plot_df["neg_log10_q"].fillna(1.0).clip(lower=0.5, upper=6.0) * 120.0

ax.scatter(
    x_vals,
    y_pos,
    s=q_sizes,
    c=point_colors,
    alpha=0.88,
    edgecolors="#1E293B",
    linewidth=1.2,
    zorder=3,
)

x_max_raw = float(np.nanmax(x_vals)) if len(x_vals) > 0 else 0.0
x_right = max(x_max_raw * 1.45, 0.015) if np.isfinite(x_max_raw) and x_max_raw > 0 else 0.02
text_offset = x_right * 0.03

for i, row in plot_df.iterrows():
    moa_brief = str(row.get("moa", "")).split(";")[0].strip()
    if len(moa_brief) > 34:
        moa_brief = moa_brief[:31].rstrip() + "..."
    phase_lbl = f"Phase {int(row['max_phase_num'])}" if row["max_phase_num"] > 0 else "Preclinical/NA"
    label_txt = f"{row['drug']}  [{phase_lbl}] — {moa_brief or 'No MOA'}"
    ax.text(
        row["abs_sim"] + text_offset,
        i,
        label_txt,
        va="center",
        fontsize=8.8,
        fontweight="bold",
        color="#1E293B",
    )

ax.set_yticks(y_pos)
ax.set_yticklabels(plot_df["bio_axis"], fontsize=9.2, fontweight="bold")
ax.invert_yaxis()
ax.set_xlim(0, x_right)
ax.set_xlabel("|Opposite Signature Similarity Score| (Bubble size ~ -log10 q-value)", fontsize=9.8, fontweight="bold")
ax.set_title(
    "Biological Axis & Clinical Stage Landscape (GSE346632: PD40 vs PD32)",
    fontsize=12,
    fontweight="bold",
    loc="left",
    pad=12,
)

legend_handles = [
    plt.Line2D([0], [0], marker="o", color="w", label="Approved (Phase 4)", markerfacecolor="#2B6CB0", markeredgecolor="#1E293B", markersize=9.5),
    plt.Line2D([0], [0], marker="o", color="w", label="Clinical / Other", markerfacecolor="#D95F02", markeredgecolor="#1E293B", markersize=9.5),
]
ax.legend(handles=legend_handles, loc="lower right", frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=8.5)

ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.grid(axis="x", linestyle="--", alpha=0.5, zorder=0)

plt.tight_layout()
plt.savefig(out_png, dpi=300)
plt.close()

print("생물학적 기전 버블 맵 저장 완료:", out_png)