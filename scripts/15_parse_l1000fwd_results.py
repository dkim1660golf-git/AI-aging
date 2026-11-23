from pathlib import Path
import re
import time
import numpy as np
import pandas as pd
import requests
from tqdm import tqdm

# 1. 입력 / 출력 경로 설정
CANDIDATE_PATHS = [
    Path("data/raw/lincs/lincs_search_results.csv"),
    Path("data/raw/lincs_search_results.csv"),
    Path("lincs_search_results.csv"),
]
in_path = next((p for p in CANDIDATE_PATHS if p.exists()), None)
if in_path is None:
    raise FileNotFoundError(
        "lincs_search_results.csv 파일을 찾을 수 없습니다.\n"
        "파일을 'data/raw/lincs/lincs_search_results.csv' 경로에 넣고 다시 실행해 주세요."
    )

out_dir = Path("results/l1000fwd")
out_dir.mkdir(parents=True, exist_ok=True)

# 2. CSV 파일 읽기 (구분자 자동 감지 포함)
df_raw = pd.read_csv(in_path, encoding="utf-8-sig", sep=None, engine="python")
print(f"[1] 원본 파일 로드 완료: {in_path}")
print(f"    - 행 수: {len(df_raw)}, 열 수: {len(df_raw.columns)}")
print(f"    - 원본 컬럼명: {list(df_raw.columns)}")

# 3. 컬럼명 표준화 (소문자 + 공백/하이픈 -> 언더스코어)
df = df_raw.copy()
df.columns = [re.sub(r"[\s\-]+", "_", str(c).strip().lower()) for c in df.columns]

def pick_col(columns, candidates):
    for cand in candidates:
        if cand in columns:
            return cand
    return None

col_sig = pick_col(df.columns, ["sig_id", "signature_id", "id"])
col_drug = pick_col(df.columns, ["drug", "pert_iname", "compound", "name", "pert_id"])
col_sim = pick_col(df.columns, ["similarity_score", "score", "scores", "similarity"])
col_pval = pick_col(df.columns, ["p_value", "pval", "pvals", "pvalue"])
col_qval = pick_col(df.columns, ["q_value", "qval", "qvals", "qvalue", "adj_p_value", "fdr"])
col_z = pick_col(df.columns, ["z_score", "zscore", "zscores"])
col_comb = pick_col(df.columns, ["combined_score", "combined_scores"])
col_dir = pick_col(df.columns, ["type", "direction", "signature_type", "cluster"])

if col_drug is None:
    raise KeyError(f"약물명 컬럼을 찾지 못했습니다. 현재 컬럼명: {list(df.columns)}")

# 표준 컬럼으로 정리
parsed = pd.DataFrame()
parsed["sig_id"] = df[col_sig].astype(str) if col_sig else ""
parsed["drug"] = df[col_drug].astype(str).str.strip()
parsed["similarity_score"] = pd.to_numeric(df[col_sim], errors="coerce") if col_sim else np.nan
parsed["p_value"] = pd.to_numeric(df[col_pval], errors="coerce") if col_pval else np.nan
parsed["q_value"] = pd.to_numeric(df[col_qval], errors="coerce") if col_qval else np.nan
parsed["z_score"] = pd.to_numeric(df[col_z], errors="coerce") if col_z else np.nan
parsed["combined_score"] = pd.to_numeric(df[col_comb], errors="coerce") if col_comb else np.nan

# Opposite(반전) 시그니처만 필터링 (CSV에 Similar/Opposite가 같이 들어있는 경우 대비)
before_n = len(parsed)
if col_dir and df[col_dir].astype(str).str.lower().str.contains("opposite").any():
    mask = df[col_dir].astype(str).str.lower().str.contains("opposite")
    parsed = parsed[mask].copy()
    print(f"[2] '{col_dir}' 컬럼 기준 Opposite 시그니처 필터링: {before_n} -> {len(parsed)}개")
elif parsed["similarity_score"].notna().any() and (parsed["similarity_score"] < 0).any():
    parsed = parsed[parsed["similarity_score"] < 0].copy()
    print(f"[2] similarity_score < 0 기준 Opposite 시그니처 필터링: {before_n} -> {len(parsed)}개")
else:
    print(f"[2] 전체 {before_n}개 행을 그대로 유지합니다.")

# 4. sig_id에서 세포주(cell_line), 시간(timepoint), 용량(dose) 추출
# 예: CPC006_A375_6H:BRD-K41859756-001-01-8:10
def parse_sig_id(sig):
    cell_line, timepoint, dose = "", "", ""
    parts = str(sig).split(":")
    if len(parts) >= 1:
        meta = parts[0].split("_")
        if len(meta) >= 3:
            cell_line = meta[1]
            timepoint = meta[2]
    if len(parts) >= 3:
        dose = parts[2]
    return pd.Series([cell_line, timepoint, dose])

parsed[["cell_line", "timepoint", "dose"]] = parsed["sig_id"].apply(parse_sig_id)

# BRD- / STK- 등 미명명 스크리닝 코드 여부 판별
CODE_PATTERN = re.compile(r"^(BRD-|STK-|MLS|NCGC|KUC|CID|LSM-|CHEMBL\d+|SCHEMBL|ZINC|SMR|HMS)", re.IGNORECASE)
parsed["is_code_only"] = parsed["drug"].apply(
    lambda x: bool(CODE_PATTERN.match(x)) or x.lower() in {"nan", "none", "", "-666"}
)

parsed_out = out_dir / "01_opposite_signatures_parsed.csv"
parsed.to_csv(parsed_out, index=False, encoding="utf-8-sig")

# 5. 약물 단위 집계
def uniq_join(series):
    vals = sorted({str(v).strip() for v in series if str(v).strip() not in {"", "nan", "None"}})
    return ", ".join(vals)

summary = (
    parsed.groupby(["drug", "is_code_only"], as_index=False)
    .agg(
        n_signatures=("sig_id", "count"),
        n_cell_lines=("cell_line", lambda s: len({v for v in s if v})),
        cell_lines=("cell_line", uniq_join),
        min_similarity_score=("similarity_score", "min"),
        best_q_value=("q_value", "min"),
        best_p_value=("p_value", "min"),
        mean_combined_score=("combined_score", "mean"),
    )
    .sort_values(["best_q_value", "min_similarity_score"], ascending=[True, True])
    .reset_index(drop=True)
)

# 6. ChEMBL REST API로 승인 상태(max_phase) 및 작용기전(MOA) 조회
CHEMBL_BASE = "https://www.ebi.ac.uk/chembl/api/data"
session = requests.Session()
session.headers.update({"Accept": "application/json"})

def fetch_chembl_info(drug_name: str) -> dict:
    empty = {
        "chembl_id": "",
        "chembl_pref_name": "",
        "max_phase": np.nan,
        "clinical_status": "Not found",
        "first_approval": "",
        "withdrawn_flag": False,
        "moa": "",
    }
    if not drug_name or CODE_PATTERN.match(drug_name):
        empty["clinical_status"] = "Code-only (skipped)"
        return empty

    mols = []
    # 1차: 공식명(pref_name) 정확 일치
    for param_key in ["pref_name__iexact", "molecule_synonyms__molecule_synonym__iexact"]:
        try:
            r = session.get(
                f"{CHEMBL_BASE}/molecule.json",
                params={param_key: drug_name, "limit": 10},
                timeout=12,
            )
            if r.ok:
                mols = r.json().get("molecules", [])
                if mols:
                    break
        except requests.RequestException:
            continue

    if not mols:
        return empty

    # max_phase가 가장 높은 항목 우선 선택
    def sort_key(m):
        mp = m.get("max_phase")
        mp_val = float(mp) if mp is not None else -1.0
        exact = 1 if str(m.get("pref_name") or "").lower() == drug_name.lower() else 0
        return (mp_val, exact)

    best = max(mols, key=sort_key)
    chembl_id = best.get("molecule_chembl_id", "")
    mp = best.get("max_phase")
    mp_val = float(mp) if mp is not None else np.nan

    if pd.isna(mp_val):
        status = "Preclinical / Unknown"
    elif mp_val >= 4:
        status = "Approved (Phase 4)"
    elif mp_val >= 3:
        status = "Phase 3"
    elif mp_val >= 2:
        status = "Phase 2"
    elif mp_val >= 1:
        status = "Phase 1"
    else:
        status = "Preclinical"

    # 작용기전(MOA) 조회
    moas = []
    if chembl_id:
        try:
            mr = session.get(
                f"{CHEMBL_BASE}/mechanism.json",
                params={"molecule_chembl_id": chembl_id, "limit": 20},
                timeout=12,
            )
            if mr.ok:
                for item in mr.json().get("mechanisms", []):
                    moa_txt = item.get("mechanism_of_action")
                    if moa_txt and moa_txt not in moas:
                        moas.append(moa_txt)
        except requests.RequestException:
            pass

    return {
        "chembl_id": chembl_id,
        "chembl_pref_name": best.get("pref_name") or "",
        "max_phase": mp_val,
        "clinical_status": status,
        "first_approval": best.get("first_approval") or "",
        "withdrawn_flag": bool(best.get("withdrawn_flag")),
        "moa": "; ".join(moas),
    }

print(f"\n[3] 고유 약물 {len(summary)}개에 대해 ChEMBL 승인 상태/MOA 조회 중...")
chembl_records = []
for drug_name in tqdm(summary["drug"], desc="ChEMBL 조회"):
    chembl_records.append(fetch_chembl_info(drug_name))
    time.sleep(0.1)

chembl_df = pd.DataFrame(chembl_records)
final_df = pd.concat([summary, chembl_df], axis=1)

# 보기 좋게 정렬 (임상 단계 높은 순 -> q-value 낮은 순 -> 유사도 더 음수인 순)
final_df["max_phase_sort"] = final_df["max_phase"].fillna(-1)
final_df = final_df.sort_values(
    ["max_phase_sort", "best_q_value", "min_similarity_score"],
    ascending=[False, True, True],
).drop(columns=["max_phase_sort"]).reset_index(drop=True)

summary_out = out_dir / "02_drug_summary_chembl.csv"
final_df.to_csv(summary_out, index=False, encoding="utf-8-sig")

# 승인 약물 + 임상 단계(Phase 1~4) 후보 별도 저장
clinical_df = final_df[final_df["max_phase"] >= 1].copy().reset_index(drop=True)
clinical_out = out_dir / "03_approved_and_clinical_candidates.csv"
clinical_df.to_csv(clinical_out, index=False, encoding="utf-8-sig")

approved_df = final_df[final_df["max_phase"] >= 4].copy().reset_index(drop=True)

# 7. 결과 출력
pd.set_option("display.max_columns", None)
pd.set_option("display.width", 160)
display_cols = [
    "drug", "clinical_status", "first_approval", "n_signatures",
    "n_cell_lines", "min_similarity_score", "best_q_value", "moa"
]

print(f"\n=== [결과 1] 승인 약물 (Approved / Phase 4): {len(approved_df)}개 ===")
if not approved_df.empty:
    print(approved_df[display_cols].to_string(index=False))
else:
    print("승인 약물로 매칭된 항목이 없습니다. (약물명 표기 차이나 코드명 위주일 수 있습니다)")

trial_df = final_df[(final_df["max_phase"] >= 1) & (final_df["max_phase"] < 4)]
print(f"\n=== [결과 2] 임상 단계 후보 (Phase 1~3): {len(trial_df)}개 ===")
if not trial_df.empty:
    print(trial_df[display_cols].to_string(index=False))
else:
    print("해당 항목 없음")

print("\n저장 완료:")
print("-", parsed_out)
print("-", summary_out)
print("-", clinical_out)