from pathlib import Path
import pandas as pd
from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats

count_path = Path("data/raw/geo/GSE346632/GSE346632_all_counts.txt.gz")
out_dir = Path("results")
out_dir.mkdir(parents=True, exist_ok=True)

# 1. 카운트 파일 불러오기
df = pd.read_csv(count_path, sep="\t", compression="infer", encoding="utf-8-sig")

sample_cols = [
    "PD32-REP1",
    "PD32-REP2",
    "PD32-REP3",
    "PD40-REP1",
    "PD40-REP2",
    "PD40-REP3",
]

# Geneid를 인덱스로 두고 샘플 카운트만 추출 후 정수형으로 정리
counts_gene_by_sample = df.set_index("Geneid")[sample_cols].copy()
counts_gene_by_sample = counts_gene_by_sample.apply(pd.to_numeric, errors="coerce").fillna(0)
counts_gene_by_sample = counts_gene_by_sample.round().astype(int)

# 2. 거의 발현되지 않는 유전자 필터링 (전체 샘플 합계 >= 10)
keep_genes = counts_gene_by_sample.sum(axis=1) >= 10
counts_filtered = counts_gene_by_sample.loc[keep_genes]

# pyDESeq2는 행=샘플, 열=유전자 형태를 입력으로 받음
counts_df = counts_filtered.T

# 3. 메타데이터(샘플 조건표) 생성
metadata = pd.DataFrame(
    {
        "condition": ["PD32", "PD32", "PD32", "PD40", "PD40", "PD40"],
    },
    index=sample_cols,
)

# 4. DESeq2 분석 실행 (PD40 vs PD32)
dds = DeseqDataSet(
    counts=counts_df,
    metadata=metadata,
    design="~condition",
)
dds.deseq2()

stat_res = DeseqStats(dds, contrast=["condition", "PD40", "PD32"])
stat_res.summary()

res_df = stat_res.results_df.copy()
res_df = res_df.reset_index()

# 인덱스 열 이름이 'Geneid' 또는 'index'일 수 있으므로 'gene_id'로 통일
if "gene_id" not in res_df.columns:
    if "Geneid" in res_df.columns:
        res_df = res_df.rename(columns={"Geneid": "gene_id"})
    elif "index" in res_df.columns:
        res_df = res_df.rename(columns={"index": "gene_id"})

res_df["ensembl_no_ver"] = res_df["gene_id"].astype(str).str.split(".").str[0]
res_df = res_df.sort_values("padj", na_position="last").reset_index(drop=True)

out_csv = out_dir / "deg_GSE346632_PD40_vs_PD32.csv"
res_df.to_csv(out_csv, index=False)

print("\n필터링 후 유전자 수:", counts_filtered.shape[0])
print("DEG 결과 저장 완료:", out_csv)
print(res_df[["gene_id", "ensembl_no_ver", "log2FoldChange", "pvalue", "padj"]].head(10))