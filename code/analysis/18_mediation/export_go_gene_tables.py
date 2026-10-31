"""Export gene names for the existing significant all-tested-background GO results."""

import argparse
import csv
import gzip
import hashlib
import json
from collections import defaultdict
from pathlib import Path


def read_table(path, delimiter="\t"):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter=delimiter))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def write_csv(path, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", required=True, type=Path)
    parser.add_argument("--workbook-json", type=Path, help="Optional intermediate for workbook formatting")
    args = parser.parse_args()
    run = args.source_run.resolve()
    source = run / "downstream_nominal_all_tested"
    significant = read_table(source / "GO_significant_BH05.tsv")
    key = lambda r: (r["screen"], r["ontology"], r["go_id"])
    complete = read_table(source / "GO_results_all.tsv.gz")
    expected = {key(r): r for r in complete if float(r["p_adjust"]) < .05}
    require(len({key(r) for r in significant}) == len(significant), "Duplicate GO associations")
    require({key(r): r for r in significant} == expected,
            "Significant export differs from the complete GO result table")
    require(not read_table(run / "downstream_nominal/GO_significant_BH05.tsv"),
            "Eligible-DEG background now has significant results; revise the report context")

    candidates = read_table(run / "report/mediation_hits_nominal_p05.csv", ",")
    by_gene = defaultdict(list)
    for row in read_table(source / "GO_gene_mapping.tsv"):
        by_gene[row["gene_id"]].append(row)
    candidate_mapping = defaultdict(lambda: defaultdict(set))
    for row in candidates:
        for mapped in by_gene[row["outcome_gene_id"]]:
            candidate_mapping[row["screen"]][mapped["entrez_id"]].add(
                (row["outcome_gene_id"], mapped["symbol"]))

    selected_terms = {r["go_id"] for r in significant}
    annotations = defaultdict(set)
    for row in read_table(source / "GO_annotation_snapshot.tsv.gz"):
        if row["go_id"] in selected_terms:
            annotations[(row["ontology"], row["go_id"])].add(row["entrez_id"])
    term_metadata = {r["go_id"]: r for r in read_table(source / "GO_term_snapshot.tsv")}
    labels = {
        "ptn_vasc_neuropil": "PTN vascular to neuropil",
        "ptn_vasc_neun": "PTN vascular to neuronal",
        "fgf1_neuropil_vasc": "FGF1 neuropil to vascular",
        "fgf2_neuropil_vasc": "FGF2 neuropil to vascular",
        "fgf1_neun_vasc": "FGF1 neuronal to vascular",
    }
    order = {sid: i for i, sid in enumerate(labels)}
    significant.sort(key=lambda r: (order[r["screen"]], float(r["p_adjust"]), r["go_id"]))
    terms, memberships = [], []
    for row in significant:
        ids = row["entrez_ids"].split(";")
        symbols = row["gene_symbols"].split(";")
        require(len(ids) == len(set(ids)) == len(symbols) == int(row["Count"]),
                f"Gene-count mismatch for {key(row)}")
        annotated = annotations[(row["ontology"], row["go_id"])]
        require(set(ids) == set(candidate_mapping[row["screen"]]) & annotated,
                f"Candidate/annotation intersection mismatch for {key(row)}")
        require(term_metadata[row["go_id"]]["description"] == row["description"],
                f"Term description differs from snapshot for {key(row)}")
        for entrez, symbol in zip(ids, symbols):
            matches = candidate_mapping[row["screen"]][entrez]
            require(len(matches) == 1, f"Ambiguous candidate mapping for {entrez}")
            gene_id, mapped_symbol = next(iter(matches))
            require(symbol == mapped_symbol, f"Gene symbol mismatch for {entrez}")
            memberships.append({
                "gene_symbol": symbol, "outcome_gene_id": gene_id, "entrez_id": entrez,
                "screen": row["screen"], "screen_label": labels[row["screen"]],
                "go_term": row["description"], "go_id": row["go_id"],
                "ontology": row["ontology"], "p_adjust": row["p_adjust"],
            })
        terms.append({
            "go_term": row["description"], "gene_symbols": ";".join(sorted(symbols)),
            "n_candidate_genes": int(row["Count"]), "screen": row["screen"],
            "screen_label": labels[row["screen"]], "go_id": row["go_id"],
            "ontology": row["ontology"], "candidate_annotated_genes": int(row["query_annotated"]),
            "background_term_genes": int(row["term_background"]),
            "background_annotated_genes": int(row["background_annotated"]),
            **{field: row[field] for field in ("GeneRatio", "BgRatio", "fold_enrichment",
                "pvalue", "p_adjust", "p_adjust_all_eligible_terms", "p_adjust_all_screens_ontologies")},
        })
    memberships.sort(key=lambda r: (r["gene_symbol"], order[r["screen"]], r["go_id"]))

    ## validate the existing FGF overlap references independently of the plotted counts.
    overlap_path = run / "downstream_nominal/gene_overlap_intersections.tsv"
    overlap = read_table(overlap_path)
    matrix_path = run / "downstream_nominal/gene_overlap_membership.tsv"
    matrix = read_table(matrix_path)
    fgf_shared = {r["gene_symbol"] for r in matrix
                  if r["fgf1_neuropil_vasc"] == r["fgf2_neuropil_vasc"] == "1"}
    two_only = next(r for r in overlap
                    if r["screens"] == "fgf1_neuropil_vasc;fgf2_neuropil_vasc")
    three = next(r for r in overlap
                 if r["screens"] == "fgf1_neuropil_vasc;fgf2_neuropil_vasc;fgf1_neun_vasc")
    require(fgf_shared == set(two_only["gene_symbols"].split(";"))
            | set(three["gene_symbols"].split(";")), "FGF overlap tables disagree")
    require(len(fgf_shared) == int(two_only["n_genes"]) + int(three["n_genes"]),
            "Exclusive FGF overlaps do not sum to the inclusive overlap")

    out = run / "report"
    term_file = out / "GO_term_genes_all_tested_BH05.csv"
    member_file = out / "GO_gene_term_membership_all_tested_BH05.csv"
    write_csv(term_file, terms)
    write_csv(member_file, memberships)
    unique_genes = sorted({r["gene_symbol"] for r in memberships})
    n_terms = len({r["go_id"] for r in terms})
    text = [
        "# Genes in significant GO associations", "",
        f"{len(terms)} screen–term associations, {n_terms} distinct GO terms and "
        f"{len(unique_genes)} distinct candidate genes. Each row lists the candidate outcome genes "
        "annotated to that term, rather than all genes annotated to the term in the background.", "",
        "Selection: the existing **all-tested target-gene background** results with "
        "BH-adjusted p < .05 within each screen/ontology. Statistics are unchanged. "
        "The eligible-DEG background has no significant terms. BP = biological process; "
        "MF = molecular function.", "",
        f"[Association CSV]({term_file.name}) · "
        f"[One gene per screen–term row]({member_file.name})", "",
        "| Screen | GO term | GO ID | Candidate genes | Count | BH-adjusted p |",
        "| --- | --- | --- | --- | ---: | ---: |",
    ]
    for row in terms:
        genes = ", ".join(f"**{g}**" for g in row["gene_symbols"].split(";"))
        text.append(f"| {row['screen_label']} | {row['go_term']} ({row['ontology']}) | "
                    f"{row['go_id']} | {genes} | {row['n_candidate_genes']} | "
                    f"{float(row['p_adjust']):.5f} |")
    text += ["", "## Interpretation", "",
        "Seven FGF2 metal-response terms contain the same **MT1G, MT1M and MT3** genes. "
        "The growth terms overlap those genes and additional candidates. "
        "Repeated annotations do not represent independent pathway evidence. "
        "These results neither establish metal exposure nor validate causal mediation.", "",
        "## FGF1/FGF2 genes in the UpSet plot", "",
        "The existing [exact-intersection table](../downstream_nominal/gene_overlap_intersections.tsv) "
        "lists gene symbols for each plotted column. The "
        "[membership matrix](../downstream_nominal/gene_overlap_membership.tsv) "
        "provides one gene per row and one membership indicator per screen.", "",
        f"The neuropil FGF1 and FGF2 screens share {len(fgf_shared)} genes: "
        f"{two_only['n_genes']} in those two screens only, plus {three['n_genes']} "
        "also in the neuronal FGF1 screen (**FAM107A, SLC24A2**). "
        "The UpSet bars show exact intersections, so these categories are separate.", "",
        "## Sources and reproduction", "",
        "Source: [original significant GO table](../downstream_nominal_all_tested/GO_significant_BH05.tsv). "
        "[Methods](../../../../code/analysis/18_mediation/DOWNSTREAM.md). "
        "CSV p-values retain the original precision. Other BH columns are separate sensitivity corrections.", "",
        "From the repository root:", "", "```bash",
        "python3 code/analysis/18_mediation/export_go_gene_tables.py \\",
        "  --source-run processed-data/18_mediation/run_2026-10-08_limma", "```", "",
    ]
    markdown = out / "GO_term_genes_all_tested_BH05.md"
    markdown.write_text("\n".join(text))
    inputs = [source / name for name in ("GO_significant_BH05.tsv", "GO_results_all.tsv.gz",
        "GO_gene_mapping.tsv", "GO_annotation_snapshot.tsv.gz", "GO_term_snapshot.tsv")]
    inputs += [run / "report/mediation_hits_nominal_p05.csv", overlap_path, matrix_path,
               run / "downstream_nominal/GO_significant_BH05.tsv"]
    fingerprint = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
    provenance = {
        "background": "all_tested", "selection": "p_adjust < 0.05 within each screen/ontology",
        "statistics": "unchanged source results", "associations": len(terms),
        "distinct_go_terms": n_terms, "distinct_candidate_genes": unique_genes,
        "gene_screen_term_rows": len(memberships),
        "fgf_neuropil_shared_genes": sorted(fgf_shared),
        "inputs_sha256": {str(p.relative_to(run)): fingerprint(p) for p in inputs},
        "outputs_sha256": {p.name: fingerprint(p) for p in (term_file, member_file, markdown)},
    }
    (out / "GO_gene_tables_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if args.workbook_json:
        args.workbook_json.write_text(json.dumps({"terms": terms, "memberships": memberships}))
    print(json.dumps({k: provenance[k] for k in ("associations", "distinct_go_terms",
        "distinct_candidate_genes", "gene_screen_term_rows")}))


if __name__ == "__main__":
    main()
