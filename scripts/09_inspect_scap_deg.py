from pathlib import Path
import pandas as pd

deg_path = Path("results/deg_GSE346632_PD40_vs_PD32_annotated.csv")
scap_path = Path("data/processed/signatures/scap_signature_starter.csv")
out_csv = Path("results/scap_starter_deg_GSE346632.csv")

deg_df = pd.read_csv(deg_path)
scap_df = pd.read_csv(scap_path)

scap_genes = (
    scap_df["gene_symbol"]
    .dropna()
    .astype(str)
    .str.strip()
    .str.upper()
    .unique()
    .tolist()
)

sub = deg_df[deg_df["symbol"].astype(str).str.strip().str.upper().isin(scap_genes)].copy()
sub["symbol_clean"] = sub["symbol"].astype(str).str.strip().str.upper()

# padj 기준 오름차순 정렬 후 symbol별 대표 행 정리
sub = sub.sort_values(["symbol_clean", "padj"], ascending=[True, True], na_position="last")

def classify_deg(row):
    if pd.isna(row["padj"]) or row["padj"] >= 0.05:
        return "NS (padj>=0.05 or NA)"
    if row["log2FoldChange"] > 0:
        return "UP_sig (LFC>0)"
    if row["log2FoldChange"] < 0:
        return "DOWN_sig (LFC<0)"
    return "FLAT"

sub["deg_status"] = sub.apply(classify_deg, axis=1)

# 경로 그룹 정보도 함께 붙이기
if "pathway_group" in scap_df.columns:
    pw_map = (
        scap_df.assign(sym=scap_df["gene_symbol"].astype(str).str.strip().str.upper())
        .drop_duplicates("sym")
        .set_index("sym")["pathway_group"]
        .to_dict()
    )
    sub["pathway_group"] = sub["symbol_clean"].map(pw_map)

cols = [
    c
    for c in [
        "symbol_clean",
        "pathway_group",
        "gene_id",
        "baseMean",
        "log2FoldChange",
        "pvalue",
        "padj",
        "deg_status",
    ]
    if c in sub.columns
]

out_df = sub[cols].rename(columns={"symbol_clean": "symbol"})
out_df.to_csv(out_csv, index=False)

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 140)
print(f"저장 완료: {out_csv} (총 {len(out_df)}행)")
print(out_df.to_string(index=False))