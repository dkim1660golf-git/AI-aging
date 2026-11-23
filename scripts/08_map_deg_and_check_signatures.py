from pathlib import Path
import pandas as pd

deg_path = Path("results/deg_GSE346632_PD40_vs_PD32.csv")
sig_dir = Path("data/processed/signatures")
out_dir = Path("results")
out_dir.mkdir(parents=True, exist_ok=True)

deg_df = pd.read_csv(deg_path)
deg_df["ensembl_no_ver"] = deg_df["gene_id"].astype(str).str.split(".").str[0]

# 1. 전체 DEG 기본 요약 (PD40 vs PD32)
valid_padj = deg_df["padj"].dropna()
sig_deg = deg_df[deg_df["padj"] < 0.05].copy()
up_deg = sig_deg[sig_deg["log2FoldChange"] > 0]
down_deg = sig_deg[sig_deg["log2FoldChange"] < 0]

print("[GSE346632 PD40 vs PD32 DEG 요약]")
print(f"- 전체 테스트 유전자 수: {len(deg_df)}")
print(f"- padj 계산된 유전자 수: {len(valid_padj)}")
print(f"- 유의미한 DEG (padj < 0.05): {len(sig_deg)}개 (상향: {len(up_deg)}개, 하향: {len(down_deg)}개)")

# 2. Ensembl ID -> Gene Symbol 매핑 시도 (mygene 패키지 활용)
mapped_df = deg_df.copy()
mapped_df["symbol"] = pd.NA

try:
    import mygene

    mg = mygene.MyGeneInfo()
    unique_ids = deg_df["ensembl_no_ver"].dropna().unique().tolist()
    res = mg.querymany(
        unique_ids,
        scopes="ensembl.gene",
        fields="symbol",
        species="human",
        returnall=False,
        verbose=False,
    )
    id_to_symbol = {}
    for item in res:
        q = item.get("query")
        sym = item.get("symbol")
        if q and sym and q not in id_to_symbol:
            id_to_symbol[q] = str(sym).strip().upper()

    mapped_df["symbol"] = mapped_df["ensembl_no_ver"].map(id_to_symbol)
    print(f"- Symbol 매핑 성공 유전자 수: {mapped_df['symbol'].notna().sum()} / {len(mapped_df)}")
except Exception as e:
    print(f"- 참고: mygene 자동 매핑은 건너뜀 ({e})")

mapped_out = out_dir / "deg_GSE346632_PD40_vs_PD32_annotated.csv"
mapped_df.to_csv(mapped_out, index=False)
print("저장 완료:", mapped_out)

# 3. Symbol이 매핑된 경우 시그니처와 교집합 확인
if mapped_df["symbol"].notna().any():
    deg_sym_set = set(mapped_df["symbol"].dropna())
    sig_up_set = set(up_deg.assign(ensembl_no_ver=up_deg["gene_id"].astype(str).str.split(".").str[0])["ensembl_no_ver"])
    # symbol 기준 유의미 상향/하향 확인
    sig_mapped = mapped_df[mapped_df["padj"] < 0.05].dropna(subset=["symbol"])
    up_sym_set = set(sig_mapped[sig_mapped["log2FoldChange"] > 0]["symbol"])
    down_sym_set = set(sig_mapped[sig_mapped["log2FoldChange"] < 0]["symbol"])

    for name, fname in [
        ("SCAP_starter", "scap_signature_starter.csv"),
        ("SenMayo", "senmayo_genes.csv"),
        ("CellAge", "cellage_genes.csv"),
        ("SASP_Atlas", "sasp_atlas_irr_fibroblast_day10.csv"),
    ]:
        fpath = sig_dir / fname
        if fpath.exists():
            sdf = pd.read_csv(fpath)
            if "gene_symbol" in sdf.columns:
                sset = set(sdf["gene_symbol"].dropna().astype(str).str.strip().str.upper())
                found_all = sorted(deg_sym_set & sset)
                found_up = sorted(up_sym_set & sset)
                found_down = sorted(down_sym_set & sset)
                print(
                    f"- {name}: 전체 매칭 {len(found_all)}/{len(sset)}개 "
                    f"| 유의 상향(padj<0.05, LFC>0): {len(found_up)}개 ({found_up[:8]}) "
                    f"| 유의 하향(padj<0.05, LFC<0): {len(found_down)}개 ({found_down[:8]})"
                )