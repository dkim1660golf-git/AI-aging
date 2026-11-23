from pathlib import Path
import numpy as np
import pandas as pd

cand_csv = Path("results/l1000fwd/04_final_repurposing_candidates.csv")
summary_csv = Path("results/l1000fwd/02_drug_summary_chembl.csv")
out_md = Path("results/l1000fwd/08_repurposing_candidates_table.md")

src_csv = cand_csv if cand_csv.exists() else summary_csv
if not src_csv.exists():
    raise FileNotFoundError("후보 또는 요약 CSV 파일을 찾지 못했습니다.")

df = pd.read_csv(src_csv)
if df.empty:
    raise ValueError("CSV 데이터가 비어 있습니다.")

df["max_phase_num"] = pd.to_numeric(df.get("max_phase", np.nan), errors="coerce")
df["min_similarity_score_num"] = pd.to_numeric(df.get("min_similarity_score", np.nan), errors="coerce")
df["best_q_value_num"] = pd.to_numeric(df.get("best_q_value", np.nan), errors="coerce")

# 실제 점수나 q-value가 있는 항목을 우선 정렬
df = df.sort_values(
    ["best_q_value_num", "min_similarity_score_num", "max_phase_num", "drug"],
    ascending=[True, True, False, True],
    na_position="last",
).reset_index(drop=True)

top_df = df.head(15).copy()

rows = []
for idx, r in top_df.iterrows():
    drug = str(r.get("drug", "")).strip()
    status = str(r.get("clinical_status", "")).strip() or "Unknown"
    phase = f"Phase {int(r['max_phase_num'])}" if pd.notna(r["max_phase_num"]) else "-"
    sim_str = f"{r['min_similarity_score_num']:.4f}" if pd.notna(r["min_similarity_score_num"]) else "-"
    q_str = f"{r['best_q_value_num']:.4g}" if pd.notna(r["best_q_value_num"]) else "-"
    moa_raw = str(r.get("moa", "")).strip()
    moa = "Not recorded" if (not moa_raw or moa_raw.lower() == "nan") else moa_raw

    rows.append(
        {
            "Rank": idx + 1,
            "Drug": drug,
            "Clinical Status": status,
            "Phase": phase,
            "Similarity": sim_str,
            "Best q-value": q_str,
            "MOA": moa,
        }
    )

table_df = pd.DataFrame(rows)

pd.set_option("display.max_columns", None)
pd.set_option("display.width", 180)
print(f"=== L1000FWD 후보 요약 테이블 (상위 {len(table_df)}개, 출처: {src_csv.name}) ===")
print(table_df.to_string(index=False))

md_lines = [
    "# L1000FWD Opposite Signature Candidate Summary",
    "",
    f"- Source file: `{src_csv}`",
    f"- Displayed rows: Top {len(table_df)}",
    "",
    "| Rank | Drug | Clinical Status | Phase | Similarity | Best q-value | MOA |",
    "| :---: | :--- | :--- | :---: | :---: | :---: | :--- |",
]
for _, r in table_df.iterrows():
    moa_escaped = r["MOA"].replace("|", "\\|")
    md_lines.append(
        f"| {r['Rank']} | {r['Drug']} | {r['Clinical Status']} | {r['Phase']} | {r['Similarity']} | {r['Best q-value']} | {moa_escaped} |"
    )

out_md.write_text("\n".join(md_lines) + "\n", encoding="utf-8")
print(f"\nMarkdown 테이블 저장 완료: {out_md}")