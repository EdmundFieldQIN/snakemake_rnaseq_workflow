#!/usr/bin/env python3

import argparse
import logging
import math
import re
from pathlib import Path
from typing import Dict, List

import numpy as np
import pandas as pd


logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge abundance quant files into counts/FPKM/TPM/TMM matrices."
    )
    parser.add_argument(
        "--est_method",
        required=True,
        choices=["featureCounts", "RSEM", "eXpress", "kallisto", "salmon"],
        help="quantification method format",
    )
    parser.add_argument(
        "--quant_files",
        default="",
        help="file containing quant file paths, one per line",
    )
    parser.add_argument(
        "--name_sample_by_basedir",
        action="store_true",
        help="use directory component as sample name",
    )
    parser.add_argument(
        "--basedir_index",
        type=int,
        default=-2,
        help="0-based index; negative allowed, default -2",
    )
    parser.add_argument(
        "--out_prefix",
        default="matrix",
        help="output prefix",
    )
    parser.add_argument(
        "files",
        nargs="*",
        help="quant files",
    )
    return parser.parse_args()


def get_colspec(est_method: str) -> Dict[str, str]:
    m = est_method.lower()
    if m == "rsem":
        return {"acc": "transcript_id", "count": "expected_count", "fpkm": "FPKM", "tpm": "TPM"}
    if m == "express":
        return {"acc": "target_id", "count": "eff_counts", "fpkm": "fpkm", "tpm": "tpm"}
    if m == "kallisto":
        return {"acc": "target_id", "count": "est_counts", "fpkm": "tpm", "tpm": "tpm"}
    if m == "salmon":
        return {"acc": "Name", "count": "NumReads", "fpkm": "TPM", "tpm": "TPM"}
    if m == "featurecounts":
        return {"acc": "gene_id", "count": "counts", "fpkm": "fpkm", "tpm": "tpm"}
    raise ValueError("Unsupported est_method")


def resolve_files(args: argparse.Namespace) -> List[str]:
    if args.quant_files:
        lines = Path(args.quant_files).read_text(encoding="utf-8").splitlines()
        files = [x.strip() for x in lines if x.strip()]
    else:
        files = list(args.files)

    if not files:
        raise ValueError("No quant files provided")
    return files


def perl_index_pick(parts: List[str], idx: int) -> str:
    pos = len(parts) + idx if idx < 0 else idx
    if pos < 0 or pos >= len(parts):
        raise ValueError(f"basedir_index out of range for path: {'/'.join(parts)}")
    return parts[pos]


def sample_name(path: str, by_basedir: bool, basedir_index: int) -> str:
    p = Path(path)
    if by_basedir:
        parts = path.split("/")
        nm = perl_index_pick(parts, basedir_index)
    else:
        nm = p.name

    nm = re.sub(r"\.(genes|isoforms)\.results$", "", nm)
    nm = re.sub(r"\.count$", "", nm)
    return nm


def load_quant(file_path: str, colspec: Dict[str, str]) -> pd.DataFrame:
    df = pd.read_csv(file_path, sep="\t", header=0, dtype=str)
    need = [colspec["acc"], colspec["count"], colspec["fpkm"], colspec["tpm"]]
    missing = [c for c in need if c not in df.columns]
    if missing:
        raise ValueError(f"{file_path} missing columns: {', '.join(missing)}")

    out = pd.DataFrame(
        {
            "feature_id": df[colspec["acc"]].astype(str),
            "count": pd.to_numeric(df[colspec["count"]], errors="coerce"),
            "fpkm": pd.to_numeric(df[colspec["fpkm"]], errors="coerce"),
            "tpm": pd.to_numeric(df[colspec["tpm"]], errors="coerce"),
        }
    )
    return out


def merge_metric(quant_list: List[pd.DataFrame], sample_names: List[str], metric: str) -> pd.DataFrame:
    merged = None
    for qdf, sname in zip(quant_list, sample_names):
        part = qdf[["feature_id", metric]].rename(columns={metric: sname})
        if merged is None:
            merged = part
        else:
            merged = merged.merge(part, on="feature_id", how="outer", sort=False)
    return merged


def calc_tmm_factors(counts: pd.DataFrame, logratio_trim: float = 0.30, sum_trim: float = 0.05) -> pd.Series:
    mat = counts.to_numpy(dtype=float)
    mat = np.nan_to_num(mat, nan=0.0)
    lib_sizes = mat.sum(axis=0)

    if mat.shape[1] < 2:
        return pd.Series(np.ones(mat.shape[1]), index=counts.columns, dtype=float)

    ref_idx = int(np.argsort(lib_sizes)[len(lib_sizes) // 2])
    ref_counts = mat[:, ref_idx]
    ref_lib = lib_sizes[ref_idx]

    factors = np.ones(mat.shape[1], dtype=float)

    for j in range(mat.shape[1]):
        if j == ref_idx:
            continue

        y = mat[:, j]
        lib_j = lib_sizes[j]

        keep = (y > 0) & (ref_counts > 0) & (lib_j > 0) & (ref_lib > 0)
        if keep.sum() < 10:
            factors[j] = 1.0
            continue

        yk = y[keep]
        rk = ref_counts[keep]

        m_vals = np.log2((yk / lib_j) / (rk / ref_lib))
        a_vals = 0.5 * np.log2((yk / lib_j) * (rk / ref_lib))

        w = 1.0 / (
            ((lib_j - yk) / (lib_j * yk)) +
            ((ref_lib - rk) / (ref_lib * rk))
        )

        m_lo, m_hi = np.quantile(m_vals, [logratio_trim, 1.0 - logratio_trim])
        a_lo, a_hi = np.quantile(a_vals, [sum_trim, 1.0 - sum_trim])

        use = (m_vals >= m_lo) & (m_vals <= m_hi) & (a_vals >= a_lo) & (a_vals <= a_hi)
        if use.sum() == 0:
            factors[j] = 1.0
            continue

        mean_m = np.sum(w[use] * m_vals[use]) / np.sum(w[use])
        factors[j] = 2.0 ** mean_m

    gm = math.exp(np.mean(np.log(factors)))
    factors = factors / gm
    return pd.Series(factors, index=counts.columns, dtype=float)


def tmm_cpm(matrix_df: pd.DataFrame) -> pd.DataFrame:
    numeric = matrix_df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    lib_sizes = numeric.sum(axis=0)
    factors = calc_tmm_factors(numeric)

    eff_lib = lib_sizes * factors
    eff_lib = eff_lib.replace(0, np.nan)

    cpm = numeric.divide(eff_lib, axis=1) * 1_000_000
    cpm = cpm.fillna(0.0)
    return cpm


def main() -> None:
    args = parse_args()
    files = resolve_files(args)
    colspec = get_colspec(args.est_method)

    sample_names = [
        sample_name(f, args.name_sample_by_basedir, args.basedir_index)
        for f in files
    ]
    if len(set(sample_names)) != len(sample_names):
        dup = sorted({x for x in sample_names if sample_names.count(x) > 1})
        raise ValueError(f"Duplicated sample names: {', '.join(dup)}")

    quant_list = [load_quant(f, colspec) for f in files]

    counts_mat = merge_metric(quant_list, sample_names, "count")
    fpkm_mat = merge_metric(quant_list, sample_names, "fpkm")
    tpm_mat = merge_metric(quant_list, sample_names, "tpm")

    counts_file = f"{args.out_prefix}.counts.matrix"
    fpkm_file = f"{args.out_prefix}.FPKM.EXPR.matrix"
    tpm_file = f"{args.out_prefix}.TPM.EXPR.matrix"
    tmm_file = f"{args.out_prefix}.TMM.EXPR.matrix"
    tpm_tmm_file = f"{args.out_prefix}.TPM.TMM.EXPR.matrix"

    counts_mat.to_csv(counts_file, sep="\t", index=False, na_rep="NA")
    fpkm_mat.to_csv(fpkm_file, sep="\t", index=False, na_rep="NA")
    tpm_mat.to_csv(tpm_file, sep="\t", index=False, na_rep="NA")

    count_numeric = counts_mat.set_index("feature_id").apply(pd.to_numeric, errors="coerce").fillna(0.0).round()
    tmm_count = tmm_cpm(count_numeric)
    tmm_out = tmm_count.reset_index()
    tmm_out.to_csv(tmm_file, sep="\t", index=False, na_rep="NA")

    tpm_numeric = tpm_mat.set_index("feature_id").apply(pd.to_numeric, errors="coerce").fillna(0.0)
    tmm_tpm = tmm_cpm(tpm_numeric)
    tmm_tpm_out = tmm_tpm.reset_index()
    tmm_tpm_out.to_csv(tpm_tmm_file, sep="\t", index=False, na_rep="NA")


if __name__ == "__main__":
    main()