from pathlib import Path
import pandas as pd

count_path = Path("data/raw/geo/GSE346632/GSE346632_all_counts.txt.gz")

df = pd.read_csv(count_path, sep="\t", compression="infer", encoding="utf-8-sig")

print("전체 행/열 수:", df.shape)
print("열 목록:", list(df.columns))
print(df.head(5))