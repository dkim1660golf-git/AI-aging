from pathlib import Path
import pandas as pd

xlsx_path = Path("data/raw/gene_sets/SASP_IrradiationFibroblast_day10_secretome_candidates_v2.xlsx")
sheet_name = "irr_day10_secretome_candidates_"
out_dir = Path("data/processed/signatures")
out_dir.mkdir(parents=True, exist_ok=True)

df = pd.read_excel(xlsx_path, sheet_name=sheet_name)

# 'Genes' 열의 세미콜론(;) 구분 유전자들을 행 단위로 펼친 뒤 정리
exploded = df.assign(
    gene_symbol=df["Genes"].fillna("").astype(str).str.split(";")
).explode("gene_symbol")

exploded["gene_symbol"] = (
    exploded["gene_symbol"]
    .astype(str)
    .str.strip()
    .str.upper()
)

clean_df = exploded[
    exploded["gene_symbol"].ne("") & exploded["gene_symbol"].ne("NAN")
].copy()

# 유전자 기호 기준 중복 제거 (필요시 Log2 Ratio 절댓값 큰 순으로 정렬 후 중복 제거 가능)
clean_df = clean_df.drop_duplicates(subset=["gene_symbol"]).reset_index(drop=True)

out_csv = out_dir / "sasp_atlas_irr_fibroblast_day10.csv"
clean_df.to_csv(out_csv, index=False)

print("저장 완료:", out_csv)
print("원본 단백질 행 수:", len(df))
print("정리된 고유 유전자 수:", len(clean_df))
print(clean_df[["gene_symbol", "Log2 Ratio", "Qvalue"]].head(10))