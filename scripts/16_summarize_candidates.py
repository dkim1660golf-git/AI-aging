from pathlib import Path
import pandas as pd

cand_path = Path("results/l1000fwd/03_approved_and_clinical_candidates.csv")
all_chembl_path = Path("results/l1000fwd/02_drug_summary_chembl.csv")
scap_deg_path = Path("results/scap_starter_deg_GSE346632.csv")
out_csv = Path("results/l1000fwd/04_final_repurposing_candidates.csv")

if not cand_path.exists():
    raise FileNotFoundError(f"후보 파일이 없습니다: {cand_path}")

cand_df = pd.read_csv(cand_path)
all_df = pd.read_csv(all_chembl_path) if all_chembl_path.exists() else cand_df.copy()

# SCAP_starter 중 PD40에서 유의 상향(padj<0.05, LFC>0)된 유전자 목록 요약
scap_up_text = ""
if scap_deg_path.exists():
    sdf = pd.read_csv(scap_deg_path)
    up_sub = sdf[(sdf["padj"] < 0.05) & (sdf["log2FoldChange"] > 0)].sort_values(
        "log2FoldChange", ascending=False
    )
    scap_up_text = ", ".join(
        [f"{r['symbol']} (LFC={r['log2FoldChange']:+.2f})" for _, r in up_sub.iterrows()]
    )

# 승인/임상 후보 우선, 없으면 전체 상위 15개 fallback
if len(cand_df) > 0:
    target_df = cand_df.copy()
    mode_label = "승인/임상(Phase 1~4) 후보"
else:
    target_df = all_df.head(15).copy()
    mode_label = "전체 상위 15개 화합물 (승인/임상 매칭 없음)"

# 우선순위 점수 계산 (임상단계 가중치 + q-value 유의성 + 반전 강도)
phase_num = pd.to_numeric(target_df.get("max_phase", 0), errors="coerce").fillna(0)
q_val = pd.to_numeric(target_df.get("best_q_value", 1), errors="coerce").fillna(1).clip(lower=1e-12)
sim_val = pd.to_numeric(target_df.get("min_similarity_score", 0), errors="coerce").fillna(0)

target_df["priority_score"] = (
    phase_num * 10.0
    + (-pd.Series(q_val).apply(lambda x: __import__("math").log10(x))).clip(upper=10)
    + sim_val.abs() * 100.0
).round(3)

target_df = target_df.sort_values("priority_score", ascending=False).reset_index(drop=True)
target_df.insert(0, "rank", target_df.index + 1)
target_df["scap_up_targets_context"] = scap_up_text

target_df.to_csv(out_csv, index=False, encoding="utf-8-sig")

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 160)
show_cols = [
    c
    for c in [
        "rank",
        "drug",
        "clinical_status",
        "max_phase",
        "n_signatures",
        "min_similarity_score",
        "best_q_value",
        "priority_score",
        "moa",
    ]
    if c in target_df.columns
]

print(f"=== [{mode_label} 요약: 총 {len(target_df)}개] ===")
print(target_df[show_cols].to_string(index=False))
if scap_up_text:
    print(f"\n[참고] GSE346632 SCAP_starter 유의 상향 유전자: {scap_up_text}")
print(f"\n저장 완료: {out_csv}")