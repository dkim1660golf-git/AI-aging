from pathlib import Path
import numpy as np
import pandas as pd

count_path = Path("data/raw/geo/GSE346632/GSE346632_all_counts.txt.gz")
deg_annot_path = Path("results/deg_GSE346632_PD40_vs_PD32_annotated.csv")
sig_dir = Path("data/processed/signatures")
out_dir = Path("results")
out_dir.mkdir(parents=True, exist_ok=True)

sample_cols = [
    "PD32-REP1",
    "PD32-REP2",
    "PD32-REP3",
    "PD40-REP1",
    "PD40-REP2",
    "PD40-REP3",
]

# 1. 매핑된 DEG 테이블에서 ensembl_no_ver -> symbol 매핑표 만들기
deg_df = pd.read_csv(deg_annot_path)
deg_df["ensembl_no_ver"] = deg_df["gene_id"].astype(str).str.split(".").str[0]

valid_map = deg_df.dropna(subset=["symbol"]).copy()
valid_map["symbol"] = valid_map["symbol"].astype(str).str.strip().str.upper()
ens_to_sym = (
    valid_map.drop_duplicates("ensembl_no_ver")
    .set_index("ensembl_no_ver")["symbol"]
    .to_dict()
)

# 2. 원본 카운트를 불러와 symbol 기준으로 정리
raw_df = pd.read_csv(count_path, sep="\t", compression="infer", encoding="utf-8-sig")
raw_df["ensembl_no_ver"] = raw_df["Geneid"].astype(str).str.split(".").str[0]
raw_df["symbol"] = raw_df["ensembl_no_ver"].map(ens_to_sym)

expr_df = raw_df.dropna(subset=["symbol"]).copy()
for c in sample_cols:
    expr_df[c] = pd.to_numeric(expr_df[c], errors="coerce").fillna(0.0)

# 같은 symbol에 여러 Ensembl 행이 있을 경우 샘플별 합계로 통합
gene_expr = expr_df.groupby("symbol", as_index=True)[sample_cols].sum()

# log2(count + 1) 변환 및 샘플 간 유전자별 z-score 계산
log_expr = np.log2(gene_expr + 1.0)
gene_mean = log_expr.mean(axis=1)
gene_std = log_expr.std(axis=1, ddof=0).replace(0, np.nan)
z_expr = log_expr.sub(gene_mean, axis=0).div(gene_std, axis=0).fillna(0.0)

# 3. 시그니처 목록 불러오기
signatures = {
    "SCAP_starter": sig_dir / "scap_signature_starter.csv",
    "SenMayo": sig_dir / "senmayo_genes.csv",
    "CellAge": sig_dir / "cellage_genes.csv",
    "SASP_Atlas": sig_dir / "sasp_atlas_irr_fibroblast_day10.csv",
}

sample_score_rows = []
deg_summary_rows = []

# symbol 기준 DEG 대표값 (padj가 가장 작은 행 기준)
deg_by_sym = (
    valid_map.sort_values(["symbol", "padj"], ascending=[True, True], na_position="last")
    .drop_duplicates("symbol")
    .set_index("symbol")
)

for sig_name, sig_path in signatures.items():
    if not sig_path.exists():
        print(f"경고: 시그니처 파일이 없어 건너뜀 -> {sig_path}")
        continue

    sdf = pd.read_csv(sig_path)
    if "gene_symbol" not in sdf.columns:
        print(f"경고: 'gene_symbol' 열이 없어 건너뜀 -> {sig_path}")
        continue

    sig_genes = sorted(
        sdf["gene_symbol"]
        .dropna()
        .astype(str)
        .str.strip()
        .str.upper()
        .unique()
        .tolist()
    )

    present_genes = [g for g in sig_genes if g in log_expr.index]
    if len(present_genes) == 0:
        continue

    # (A) 샘플별 log2(count+1) 평균 및 z-score 평균
    mean_log_by_sample = log_expr.loc[present_genes].mean(axis=0)
    mean_z_by_sample = z_expr.loc[present_genes].mean(axis=0)

    for s in sample_cols:
        cond = "PD32" if s.startswith("PD32") else "PD40"
        sample_score_rows.append(
            {
                "signature": sig_name,
                "sample": s,
                "condition": cond,
                "n_genes_scored": len(present_genes),
                "mean_log2p1": round(float(mean_log_by_sample[s]), 4),
                "mean_zscore": round(float(mean_z_by_sample[s]), 4),
            }
        )

    # (B) 시그니처별 DEG 요약
    sub_deg = deg_by_sym.reindex([g for g in sig_genes if g in deg_by_sym.index])
    sig_only = sub_deg[sub_deg["padj"] < 0.05]
    n_up = int((sig_only["log2FoldChange"] > 0).sum())
    n_down = int((sig_only["log2FoldChange"] < 0).sum())
    mean_lfc = float(sub_deg["log2FoldChange"].mean()) if len(sub_deg) > 0 else np.nan

    deg_summary_rows.append(
        {
            "signature": sig_name,
            "total_sig_genes": len(sig_genes),
            "matched_in_expr": len(present_genes),
            "mean_log2FC": round(mean_lfc, 4) if pd.notna(mean_lfc) else np.nan,
            "up_sig_genes": n_up,
            "down_sig_genes": n_down,
        }
    )

sample_scores_df = pd.DataFrame(sample_score_rows)
deg_summary_df = pd.DataFrame(deg_summary_rows)

sample_out = out_dir / "sample_signature_scores.csv"
summary_out = out_dir / "signature_deg_summary.csv"

sample_scores_df.to_csv(sample_out, index=False)
deg_summary_df.to_csv(summary_out, index=False)

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 140)

print("=== [1] 시그니처별 DEG 요약 ===")
print(deg_summary_df.to_string(index=False))

print("\n=== [2] 샘플별 z-score 점수 피벗 (PD32 vs PD40) ===")
pivot_z = sample_scores_df.pivot(index="signature", columns="sample", values="mean_zscore")
pivot_z["PD32_mean_z"] = pivot_z[["PD32-REP1", "PD32-REP2", "PD32-REP3"]].mean(axis=1).round(4)
pivot_z["PD40_mean_z"] = pivot_z[["PD40-REP1", "PD40-REP2", "PD40-REP3"]].mean(axis=1).round(4)
pivot_z["delta_PD40_minus_PD32"] = (pivot_z["PD40_mean_z"] - pivot_z["PD32_mean_z"]).round(4)
print(pivot_z.to_string())

print("\n저장 완료:")
print("-", sample_out)
print("-", summary_out)