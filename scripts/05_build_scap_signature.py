from pathlib import Path
import pandas as pd

sig_dir = Path("data/processed/signatures")
sig_dir.mkdir(parents=True, exist_ok=True)

# 1. SCAPs 대표 경로 시그니처 스타터 정의
scap_rows = [
    # BCL-2 패밀리 (세포자멸사 억제 생존 방패)
    {"gene_symbol": "BCL2", "pathway_group": "BCL2_family", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "BCL2L1", "pathway_group": "BCL2_family", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "BCL2L2", "pathway_group": "BCL2_family", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "MCL1", "pathway_group": "BCL2_family", "direction": "up", "source": "SCAP_curated"},
    # PI3K / AKT / mTOR 생존 신호
    {"gene_symbol": "PIK3CA", "pathway_group": "PI3K_AKT", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "AKT1", "pathway_group": "PI3K_AKT", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "MTOR", "pathway_group": "PI3K_AKT", "direction": "up", "source": "SCAP_curated"},
    # Serpin / PAI 생존 및 섬유화 관련
    {"gene_symbol": "SERPINE1", "pathway_group": "SERPIN_PAI", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "SERPINE2", "pathway_group": "SERPIN_PAI", "direction": "up", "source": "SCAP_curated"},
    # HSP90 샤페론 생존 경로
    {"gene_symbol": "HSP90AA1", "pathway_group": "HSP90", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "HSP90AB1", "pathway_group": "HSP90", "direction": "up", "source": "SCAP_curated"},
    # Ephrin / 의존성 수용체 및 생존 신호
    {"gene_symbol": "EFNB1", "pathway_group": "EPH_EFNB", "direction": "up", "source": "SCAP_curated"},
    {"gene_symbol": "EFNB2", "pathway_group": "EPH_EFNB", "direction": "up", "source": "SCAP_curated"},
]

scap_df = pd.DataFrame(scap_rows)
scap_df["gene_symbol"] = scap_df["gene_symbol"].str.strip().str.upper()
scap_df = scap_df.drop_duplicates(subset=["gene_symbol", "direction"]).reset_index(drop=True)

scap_out = sig_dir / "scap_signature_starter.csv"
scap_df.to_csv(scap_out, index=False)
print(f"저장 완료: {scap_out} ({len(scap_df)}개 유전자)")

# 2. 기존 정리된 유전자 세트와 교집합 간단 점검
sets = {"SCAP_starter": set(scap_df["gene_symbol"])}
for name, fname in [
    ("SenMayo", "senmayo_genes.csv"),
    ("CellAge", "cellage_genes.csv"),
    ("SASP_Atlas", "sasp_atlas_irr_fibroblast_day10.csv"),
]:
    fpath = sig_dir / fname
    if fpath.exists():
        df_tmp = pd.read_csv(fpath)
        if "gene_symbol" in df_tmp.columns:
            sets[name] = set(df_tmp["gene_symbol"].dropna().astype(str).str.strip().str.upper())

print("\n[유전자 세트별 고유 유전자 수]")
for k, v in sets.items():
    print(f"- {k}: {len(v)}개")

if "SCAP_starter" in sets:
    scap_set = sets["SCAP_starter"]
    for k, v in sets.items():
        if k != "SCAP_starter":
            overlap = sorted(scap_set & v)
            print(f"- SCAP_starter ∩ {k} ({len(overlap)}개): {overlap}")