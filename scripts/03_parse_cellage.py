from pathlib import Path
import pandas as pd

# 1. 파일 경로 설정 (실제 파일명으로 수정하세요)
tsv_path = Path("data/raw/gene_sets/cellage_genes.tsv")
out_dir = Path("data/processed/signatures")
out_dir.mkdir(parents=True, exist_ok=True)

# 2. TSV 파일 읽기
df = pd.read_csv(tsv_path, sep="\t", encoding="utf-8-sig")

# 3. 실제 열 이름('Gene symbol') 기준으로 유전자 기호 정리
gene_col = "Gene symbol"

clean_df = df.copy()
clean_df["gene_symbol"] = (
    clean_df[gene_col]
    .astype(str)
    .str.strip()
    .str.upper()
)

# 빈 값이나 'NAN' 문자열 제거
clean_df = clean_df[
    clean_df["gene_symbol"].ne("") & clean_df["gene_symbol"].ne("NAN")
]

# 유전자 기호 기준 중복 제거
clean_df = clean_df.drop_duplicates(subset=["gene_symbol"]).reset_index(drop=True)

# 4. CSV로 저장
out_csv = out_dir / "cellage_genes.csv"
clean_df.to_csv(out_csv, index=False)

print("저장 완료:", out_csv)
print("원본 행 수:", len(df))
print("정리된 고유 유전자 수:", len(clean_df))
print(clean_df[["gene_symbol"]].head(10))