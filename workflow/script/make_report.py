#!/usr/bin/env python3
"""
Generate an HTML summary report from the RNA-seq Snakemake pipeline output.

Sections:
  1. Sample QC quality (from fastp JSON)
  2. Alignment statistics (from HISAT2 logs)
  3. Expression matrix summary (from count matrix)
  4. Differential expression results (from DESeq2 output)
"""

import argparse
import json
import os
import sys
from datetime import datetime
from pathlib import Path


# ---------------------------------------------------------------------------
# Data collectors
# ---------------------------------------------------------------------------

def collect_fastp_qc(workdir: str) -> list[dict]:
    """Parse fastp JSON files and return per-sample QC records."""
    rundir = Path(workdir) / "01cleandata"
    records = []
    for jf in sorted(rundir.glob("*.fastp.json")):
        sample = jf.name.replace(".fastp.json", "")
        with open(jf) as fh:
            data = json.load(fh)
        s = data["summary"]
        bf = s["before_filtering"]
        af = s["after_filtering"]
        fr = data.get("filtering_result", {})
        dup = data.get("duplication", {})
        records.append({
            "sample": sample,
            "raw_reads": bf["total_reads"],
            "clean_reads": af["total_reads"],
            "retention": af["total_reads"] / max(bf["total_reads"], 1) * 100,
            "q20_before": bf.get("q20_rate", 0) * 100,
            "q30_before": bf.get("q30_rate", 0) * 100,
            "q20_after": af.get("q20_rate", 0) * 100,
            "q30_after": af.get("q30_rate", 0) * 100,
            "gc_content": af.get("gc_content", 0) * 100,
            "dup_rate": dup.get("rate", 0) * 100,
            "low_quality": fr.get("low_quality_reads", 0),
            "too_many_N": fr.get("too_many_N_reads", 0),
            "adapter_dimer": fr.get("adapter_dimer_reads", 0),
        })
    return records


def collect_alignment_stats(workdir: str) -> list[dict]:
    """Parse HISAT2 log files and return per-sample alignment records."""
    rundir = Path(workdir) / "03aligned"
    records = []

    def _pct(val: str) -> float:
        return float(val.strip().strip("()").rstrip("%"))

    for lf in sorted(rundir.glob("*.hisat2.log")):
        sample = lf.name.replace(".hisat2.log", "")
        text = lf.read_text()

        rec = {"sample": sample}
        # Overall alignment rate
        for line in text.splitlines():
            if "Overall alignment rate" in line:
                rec["overall_rate"] = _pct(line.split()[-1])
        # Concordant pairs
        for line in text.splitlines():
            if "Aligned concordantly or discordantly 0 time" in line:
                parts = line.strip().split()
                rec["unaligned_pairs_pct"] = _pct(parts[-1])
            if "Aligned concordantly 1 time" in line:
                parts = line.strip().split()
                rec["unique_pairs_pct"] = _pct(parts[-1])
            if "Aligned concordantly >1 times" in line:
                parts = line.strip().split()
                rec["multi_pairs_pct"] = _pct(parts[-1])
        records.append(rec)
    return records


def collect_expression_summary(workdir: str) -> dict:
    """Summarise the gene count matrix."""
    mat_path = Path(workdir) / "05expression" / "genes.counts.matrix"
    if not mat_path.is_file():
        return {}
    with open(mat_path) as fh:
        header = fh.readline().strip().split("\t")
    samples = header[1:]  # first column is feature_id
    # Count expressed genes per sample (count > 0)
    expressed = {s: 0 for s in samples}
    total_counts = {s: 0 for s in samples}
    n_genes = 0
    with open(mat_path) as fh:
        fh.readline()  # skip header
        for line in fh:
            n_genes += 1
            parts = line.strip().split("\t")
            for i, s in enumerate(samples, start=1):
                c = int(parts[i])
                if c > 0:
                    expressed[s] += 1
                total_counts[s] += c
    return {
        "n_genes": n_genes,
        "n_samples": len(samples),
        "samples": samples,
        "expressed": expressed,
        "total_counts": total_counts,
    }


def collect_de_summary(workdir: str) -> list[dict]:
    """Summarise DE results per contrast."""
    dedir = Path(workdir) / "06de_analysis" / "DESeq2_Result"
    records = []
    for fp in sorted(dedir.glob("*DE_results")):
        name = fp.name
        # e.g. genes.counts.matrix.4_vs_WT.DESeq2.DE_results
        contrast = name.replace("genes.counts.matrix.", "").replace(".DESeq2.DE_results", "")
        up = down = total = padj_sig = 0
        with open(fp) as fh:
            fh.readline()  # header
            for line in fh:
                parts = line.strip().split("\t")
                if len(parts) < 11:
                    continue
                try:
                    padj = float(parts[10])
                    log2fc = float(parts[6])
                except ValueError:
                    continue
                total += 1
                if padj < 0.05:
                    padj_sig += 1
                    if log2fc > 1:
                        up += 1
                    elif log2fc < -1:
                        down += 1
        records.append({
            "contrast": contrast,
            "total_genes": total,
            "padj_sig": padj_sig,
            "up": up,
            "down": down,
        })
    return records


def collect_quant_stats(workdir: str) -> list[dict]:
    """Parse featureCounts summary logs."""
    rundir = Path(workdir) / "04quantification"
    records = []
    for lf in sorted(rundir.glob("*.log")):
        sample = lf.name.replace(".log", "")
        rec = {"sample": sample}
        for line in lf.read_text().splitlines():
            if "\t" in line:
                k, v = line.strip().split("\t", 1)
                rec[k] = v
        records.append(rec)
    return records


# ---------------------------------------------------------------------------
# HTML builder
# ---------------------------------------------------------------------------

CSS = """
body { font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;
       max-width: 1200px; margin: 0 auto; padding: 20px; color: #333; background: #fafafa; }
h1 { border-bottom: 3px solid #2c3e50; padding-bottom: 10px; color: #2c3e50; }
h2 { color: #2c3e50; margin-top: 40px; border-bottom: 2px solid #bdc3c7; padding-bottom: 6px; }
h3 { color: #555; margin-top: 24px; }
table { border-collapse: collapse; width: 100%; margin: 12px 0 24px; background: #fff;
        box-shadow: 0 1px 3px rgba(0,0,0,0.08); }
th { background: #2c3e50; color: #fff; padding: 8px 12px; text-align: right; font-weight: 500; }
th:first-child { text-align: left; }
td { padding: 6px 12px; border-bottom: 1px solid #eee; text-align: right; }
td:first-child { text-align: left; font-weight: 500; }
tr:hover { background: #f5f7fa; }
.number { font-family: 'SF Mono', 'Monaco', 'Consolas', monospace; }
.summary-box { display: flex; flex-wrap: wrap; gap: 16px; margin: 16px 0; }
.summary-item { background: #fff; border-left: 4px solid #2c3e50; padding: 12px 20px;
                 box-shadow: 0 1px 3px rgba(0,0,0,0.08); min-width: 180px; }
.summary-item .label { font-size: 0.85em; color: #888; }
.summary-item .value { font-size: 1.6em; font-weight: 600; color: #2c3e50; }
.up { color: #c0392b; font-weight: 600; }
.down { color: #2980b9; font-weight: 600; }
.sig { color: #27ae60; font-weight: 600; }
.footer { margin-top: 40px; padding-top: 12px; border-top: 1px solid #ddd;
          font-size: 0.85em; color: #999; text-align: center; }
"""


def html_page(title: str, body: str) -> str:
    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<title>{title}</title>
<style>{CSS}</style>
</head>
<body>
<h1>{title}</h1>
{body}
<div class="footer">Generated {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</div>
</body>
</html>"""


def fmt_f(val: float, decimals: int = 2) -> str:
    return f"{val:.{decimals}f}"


def fmt_n(val: int) -> str:
    return f"{val:,}"


def build_qc_section(records: list[dict]) -> str:
    if not records:
        return "<h2>1. Sample QC</h2><p>No fastp data found.</p>"

    rows = []
    for r in records:
        rows.append(f"""<tr>
<td>{r['sample']}</td>
<td class="number">{fmt_n(r['raw_reads'])}</td>
<td class="number">{fmt_n(r['clean_reads'])}</td>
<td class="number">{fmt_f(r['retention'], 1)}%</td>
<td class="number">{fmt_f(r['q20_before'], 1)}%</td>
<td class="number">{fmt_f(r['q30_before'], 1)}%</td>
<td class="number">{fmt_f(r['q20_after'], 1)}%</td>
<td class="number">{fmt_f(r['q30_after'], 1)}%</td>
<td class="number">{fmt_f(r['gc_content'], 1)}%</td>
<td class="number">{fmt_f(r['dup_rate'], 2)}%</td>
</tr>""")

    return f"""<h2>1. Sample QC (fastp)</h2>
<table>
<thead><tr>
<th>Sample</th><th>Raw Reads</th><th>Clean Reads</th><th>Retention</th>
<th>Q20 Before</th><th>Q30 Before</th><th>Q20 After</th><th>Q30 After</th>
<th>GC Content</th><th>Duplication</th>
</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>"""


def build_alignment_section(records: list[dict]) -> str:
    if not records:
        return "<h2>2. Alignment Statistics</h2><p>No HISAT2 log data found.</p>"

    rows = []
    for r in records:
        overall = r.get("overall_rate", 0)
        unique = r.get("unique_pairs_pct", 0)
        multi = r.get("multi_pairs_pct", 0)
        unaligned = r.get("unaligned_pairs_pct", 0)
        rows.append(f"""<tr>
<td>{r['sample']}</td>
<td class="number">{fmt_f(overall, 2)}%</td>
<td class="number">{fmt_f(unique, 2)}%</td>
<td class="number">{fmt_f(multi, 2)}%</td>
<td class="number">{fmt_f(unaligned, 2)}%</td>
</tr>""")

    return f"""<h2>2. Alignment Statistics (HISAT2)</h2>
<table>
<thead><tr>
<th>Sample</th><th>Overall Alignment</th><th>Unique Concordant</th>
<th>Multi-mapped Concordant</th><th>Unaligned Concordant</th>
</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>"""


def build_expression_section(expr: dict) -> str:
    if not expr:
        return "<h2>3. Expression Matrix</h2><p>No count matrix found.</p>"

    n_genes = expr["n_genes"]
    n_samples = expr["n_samples"]
    rows = []
    for s in expr["samples"]:
        total = expr["total_counts"].get(s, 0)
        exp = expr["expressed"].get(s, 0)
        pct = exp / max(n_genes, 1) * 100
        rows.append(f"""<tr>
<td>{s}</td>
<td class="number">{fmt_n(total)}</td>
<td class="number">{fmt_n(exp)}</td>
<td class="number">{fmt_f(pct, 1)}%</td>
</tr>""")

    return f"""<h2>3. Expression Matrix</h2>
<div class="summary-box">
<div class="summary-item"><div class="label">Total Genes</div><div class="value">{fmt_n(n_genes)}</div></div>
<div class="summary-item"><div class="label">Samples</div><div class="value">{n_samples}</div></div>
</div>
<h3>Per-Sample Count Summary</h3>
<table>
<thead><tr><th>Sample</th><th>Total Counts</th><th>Expressed Genes (count &gt; 0)</th><th>% of All Genes</th></tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>"""


def build_de_section(records: list[dict]) -> str:
    if not records:
        return "<h2>4. Differential Expression</h2><p>No DE results found.</p>"

    rows = []
    for r in records:
        rows.append(f"""<tr>
<td>{r['contrast']}</td>
<td class="number">{fmt_n(r['total_genes'])}</td>
<td class="number sig">{fmt_n(r['padj_sig'])}</td>
<td class="number up">{fmt_n(r['up'])}</td>
<td class="number down">{fmt_n(r['down'])}</td>
</tr>""")

    return f"""<h2>4. Differential Expression (DESeq2)</h2>
<table>
<thead><tr>
<th>Contrast</th><th>Genes Tested</th><th>padj &lt; 0.05</th>
<th class="up">Up (log2FC &gt; 1)</th><th class="down">Down (log2FC &lt; -1)</th>
</tr></thead>
<tbody>{"".join(rows)}</tbody>
</table>"""


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Generate RNA-seq HTML report")
    parser.add_argument("--workdir", "-w", default=".",
                        help="Pipeline working directory (default: .)")
    parser.add_argument("--out", "-o", default="report.html",
                        help="Output HTML file path")
    args = parser.parse_args()

    wd = args.workdir
    if not os.path.isdir(wd):
        print(f"Error: workdir not found: {wd}", file=sys.stderr)
        sys.exit(1)

    # Collect data
    qc_records = collect_fastp_qc(wd)
    aln_records = collect_alignment_stats(wd)
    expr_summary = collect_expression_summary(wd)
    de_records = collect_de_summary(wd)

    # Build sections
    sections = [
        build_qc_section(qc_records),
        build_alignment_section(aln_records),
        build_expression_section(expr_summary),
        build_de_section(de_records),
    ]
    body = "\n".join(sections)
    html = html_page("RNA-seq Analysis Report", body)

    out_path = os.path.join(wd, args.out) if not os.path.isabs(args.out) else args.out
    with open(out_path, "w", encoding="utf-8") as fh:
        fh.write(html)
    print(f"Report written to {out_path}")


if __name__ == "__main__":
    main()
