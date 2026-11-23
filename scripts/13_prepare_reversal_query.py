from pathlib import Path
import pandas as pd

deg_annot_path = Path("results/deg_GSE346632_PD40_vs_PD32_annotated.csv")
scap_deg_path = Path("results/scap_starter_deg_GSE346632.csv")
out_dir = Path("results/reversal_query")
out_dir.mkdir(parents=True, exist_ok=True)

deg_df = pd.read_csv(deg_annot_path)
valid_deg = deg_df.dropna(subset=["symbol", "log2FoldChange", "padj"]).copy()
valid_deg["symbol"] = valid_deg["symbol"].astype(str).str.strip().str.upper()

# symbol 기준 가장 유의미한(padj 최소) 행 선택
deg_unique = (
    valid_deg.sort_values(["symbol", "padj", "pvalue"], ascending=[True, True, True], na_position="last")
    .drop_duplicates("symbol")
    .reset_index(drop=True)
)

# 1. 전체 DEG 기준 유의미한 UP / DOWN 유전자 (padj < 0.05)
sig_deg = deg_unique[deg_unique["padj"] < 0.05].copy()
up_all = sig_deg[sig_deg["log2FoldChange"] > 0].sort_values("log2FoldChange", ascending=False)
down_all = sig_deg[sig_deg["log2FoldChange"] < 0].sort_values("log2FoldChange", ascending=True)

# 2. 상위 강조 쿼리용 (예: |log2FC| >= 1.0 & padj < 0.05)
strong_up = up_all[up_all["log2FoldChange"] >= 1.0].copy()
strong_down = down_all[down_all["log2FoldChange"] <= -1.0].copy()

# 3. SCAP_starter 중 유의미한 UP / DOWN 목록
scap_up_syms = []
scap_down_syms = []
if scap_deg_path.exists():
    scap_df = pd.read_csv(scap_deg_path)
    scap_sig = scap_df[scap_df["padj"] < 0.05].copy()
    scap_up_syms = sorted(scap_sig[scap_sig["log2FoldChange"] > 0]["symbol"].dropna().unique().tolist())
    scap_down_syms = sorted(scap_sig[scap_sig["log2FoldChange"] < 0]["symbol"].dropna().unique().tolist())

up_all[["symbol", "gene_id", "log2FoldChange", "padj"]].to_csv(out_dir / "query_up_all_sig.csv", index=False)
down_all[["symbol", "gene_id", "log2FoldChange", "padj"]].to_csv(out_dir / "query_down_all_sig.csv", index=False)
strong_up[["symbol", "gene_id", "log2FoldChange", "padj"]].to_csv(out_dir / "query_up_strong_lfc1.csv", index=False)
strong_down[["symbol", "gene_id", "log2FoldChange", "padj"]].to_csv(out_dir / "query_down_strong_lfc1.csv", index=False)

print("[반전 스크리닝용 노화 쿼리 시그니처 요약 (PD40 vs PD32)]")
print(f"- 전체 유의 UP (padj<0.05, LFC>0): {len(up_all)}개")
print(f"- 전체 유의 DOWN (padj<0.05, LFC<0): {len(down_all)}개")
print(f"- 강한 변화 UP (padj<0.05, LFC>=+1.0): {len(strong_up)}개")
print(f"- 강한 변화 DOWN (padj<0.05, LFC<=-1.0): {len(strong_down)}개")
print(f"- SCAP_starter 유의 UP ({len(scap_up_syms)}개): {scap_up_syms}")
print(f"- SCAP_starter 유의 DOWN ({len(scap_down_syms)}개): {scap_down_syms}")
print(f"\n저장 완료 폴더: {out_dir}")