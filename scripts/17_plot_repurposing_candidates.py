from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

in_csv = Path("results/l1000fwd/04_final_repurposing_candidates.csv")
out_png = Path("results/l1000fwd/05_repurposing_candidates_summary.png")
TOP_N = 10

if not in_csv.exists():
    raise FileNotFoundError(f"입력 파일이 없습니다: {in_csv}")

df = pd.read_csv(in_csv)
if df.empty:
    raise ValueError("후보 데이터가 비어 있어 그래프를 그릴 수 없습니다.")

df["min_similarity_score"] = pd.to_numeric(df["min_similarity_score"], errors="coerce").fillna(0.0)
df["best_q_value"] = pd.to_numeric(df["best_q_value"], errors="coerce")

if "priority_score" in df.columns:
    df["priority_score_num"] = pd.to_numeric(df["priority_score"], errors="coerce")
    df = df.sort_values(
        ["priority_score_num", "best_q_value"], ascending=[False, True], na_position="last"
    ).reset_index(drop=True)

# 너무 많은 항목이 한 번에 겹치지 않도록 상위 TOP_N개만 선택
df = df.head(TOP_N).copy().reset_index(drop=True)

valid_q = df["best_q_value"].where(df["best_q_value"] > 0).clip(lower=1e-12, upper=1.0).fillna(1.0)
df["abs_sim"] = df["min_similarity_score"].abs().replace([np.inf, -np.inf], np.nan).fillna(0.0)
df["neg_log10_q"] = (-np.log10(valid_q)).replace([np.inf, -np.inf], np.nan).fillna(0.0)


def shorten_moa(val, max_len=38):
    if pd.isna(val) or str(val).strip() == "":
        return "Unknown MOA"
    first_moa = str(val).split(";")[0].strip()
    if len(first_moa) > max_len:
        return first_moa[: max_len - 3].rstrip() + "..."
    return first_moa


df["moa_short"] = df["moa"].apply(shorten_moa) if "moa" in df.columns else "Unknown MOA"
df["status_clean"] = (
    df["clinical_status"].fillna("Unknown").replace("", "Unknown")
    if "clinical_status" in df.columns
    else "Unknown"
)

plt.rcParams["font.sans-serif"] = "DejaVu Sans"
plt.rcParams["axes.unicode_minus"] = False

fig_height = max(4.4, 0.55 * len(df) + 1.5)
fig, axes = plt.subplots(1, 2, figsize=(13.0, fig_height))
fig.patch.set_facecolor("white")

palette = ["#2B6CB0" if i % 2 == 0 else "#D95F02" for i in range(len(df))]
y_pos = np.arange(len(df))
left_labels = [f"{r['drug']} ({r['status_clean']})" for _, r in df.iterrows()]

# [Panel A] Opposite signature 반전 강도
axes[0].barh(y_pos, df["abs_sim"], height=0.52, color=palette, edgecolor="none")
axes[0].set_yticks(y_pos)
axes[0].set_yticklabels(left_labels, fontsize=9.5, fontweight="bold")
axes[0].invert_yaxis()
axes[0].set_xlabel("|Similarity Score| (Opposite Signature)", fontsize=9.5, fontweight="bold")
axes[0].set_title(f"A. Reversal Magnitude (Top {len(df)})", fontsize=11, fontweight="bold", loc="left", pad=8)

raw_max_sim = float(np.nanmax(df["abs_sim"].values)) if len(df) > 0 else 0.0
max_sim = raw_max_sim if np.isfinite(raw_max_sim) and raw_max_sim > 0 else 0.03
axes[0].set_xlim(0, max_sim * 1.28)
for i, v in enumerate(df["abs_sim"]):
    axes[0].text(v + max_sim * 0.025, i, f"{v:.4f}", va="center", fontsize=8.5, fontweight="bold", color="#1E293B")

# [Panel B] 통계 유의성 (-log10 q-value) & 작용기전(MOA)
axes[1].barh(y_pos, df["neg_log10_q"], height=0.52, color=palette, edgecolor="none")
q_cutoff = -np.log10(0.05)
axes[1].axvline(q_cutoff, color="#64748B", linestyle="--", linewidth=1.2, label="q = 0.05 threshold")
axes[1].set_yticks(y_pos)
axes[1].set_yticklabels(df["moa_short"], fontsize=9, fontweight="bold")
axes[1].invert_yaxis()
axes[1].set_xlabel("-log10(best q-value)", fontsize=9.5, fontweight="bold")
axes[1].set_title("B. Significance & Primary MOA", fontsize=11, fontweight="bold", loc="left", pad=8)

raw_max_q = float(np.nanmax(df["neg_log10_q"].values)) if len(df) > 0 else 0.0
max_q = max(raw_max_q if np.isfinite(raw_max_q) else 0.0, q_cutoff * 1.15)
axes[1].set_xlim(0, max_q * 1.28)
for i, v in enumerate(df["neg_log10_q"]):
    axes[1].text(v + max_q * 0.02, i, f"{v:.2f}", va="center", fontsize=8.5, fontweight="bold", color="#1E293B")

axes[1].legend(loc="lower right", frameon=True, facecolor="white", edgecolor="#CBD5E1", fontsize=8.5)

for ax in axes:
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(axis="x", linestyle=":", alpha=0.5)

fig.suptitle(
    "L1000FWD Opposite Signature Candidates (GSE346632: PD40 vs PD32)",
    fontsize=12.5,
    fontweight="bold",
    y=0.98,
)
plt.tight_layout(rect=[0, 0, 1, 0.93], w_pad=3.0)
plt.savefig(out_png, dpi=300)
plt.close()

print(f"후보 요약 그래프 저장 완료 (상위 {len(df)}개): {out_png}")