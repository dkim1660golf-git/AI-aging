from pathlib import Path
import json
import pandas as pd

json_path = Path("data/raw/gene_sets/SAUL_SEN_MAYO.v2026.2.Hs.json")
out_dir = Path("data/processed/signatures")
out_dir.mkdir(parents=True, exist_ok=True)

with json_path.open("r", encoding="utf-8") as f:
    data = json.load(f)

entry = data["SAUL_SEN_MAYO"]
genes = [g.strip().upper() for g in entry["geneSymbols"]]
genes = list(dict.fromkeys(genes))  # 순서 유지하며 중복 제거

df = pd.DataFrame({
    "gene_symbol": genes,
    "gene_set": "SAUL_SEN_MAYO",
    "systematic_name": entry.get("systematicName", ""),
    "pmid": entry.get("pmid", ""),
    "source": "MSigDB_SenMayo"
})

out_csv = out_dir / "senmayo_genes.csv"
df.to_csv(out_csv, index=False)

print("저장 완료:", out_csv)
print("추출된 유전자 수:", len(df))
print(df.head(10))