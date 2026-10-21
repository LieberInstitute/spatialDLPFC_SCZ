"""Export nominal-p screening candidates from an unchanged mediation run."""

import argparse
import csv
import gzip
import hashlib
import json
from collections import Counter
from pathlib import Path


def read_table(path):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_table(path, rows, fields, delimiter="\t"):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter=delimiter, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", required=True, type=Path)
    args = parser.parse_args()
    run = args.source_run.resolve()
    source = run / "primary/all_pairs.tsv.gz"
    source_sha = hashlib.sha256(source.read_bytes()).hexdigest()
    rows = read_table(source)
    manifest = read_table(Path(__file__).with_name("screens.tsv"))
    screen_order = {r["screen_id"]: i for i, r in enumerate(manifest)}
    key = lambda r: (r["screen_id"], r["gene_id"])
    signatures = {r["run_signature"] for r in rows}
    assert len(signatures) == 1, "Mixed source-run signatures"
    assert len({key(r) for r in rows}) == len(rows), "Duplicate source pairs"
    assert {r["screen_id"] for r in rows} == set(screen_order), "Incomplete screen family"
    assert all(r["pooled_family_complete"] == "TRUE" for r in rows)
    for screen in manifest:
        selected = [r for r in rows if r["screen_id"] == screen["screen_id"]]
        assert {int(r["n_samples"]) for r in selected} == {int(screen["expected_samples"])}
        assert {int(r["n_donors"]) for r in selected} == {int(screen["expected_donors"])}

    ## preserve the historical and matched baseline gates, and significance-loss rule.
    ## matched mediator diagnosis p-values remain annotations, not selection gates.
    eligible = [r for r in rows if float(r["historical_p"]) < .05
                and float(r["c_p"]) < .05 and float(r["cprime_p"]) >= .05]
    previous = [r for r in eligible if float(r["b_q"]) < .10]
    nominal = [r for r in eligible if float(r["b_p"]) < .05]
    previous_keys = {key(r) for r in previous}
    assert previous_keys == {key(r) for r in rows if r["screen_hit"] == "TRUE"}, \
        "Source hit flags do not match the expected previous screening policy"
    reference = run / "report/mediation_hits_92.tab"
    if reference.exists():
        archived = read_table(reference)
        assert previous_keys == {(r["screen"], r["outcome_gene_id"]) for r in archived}, \
            "Previous pair membership differs from the archived table"
    assert previous_keys <= {key(r) for r in nominal}, "Nominal filter unexpectedly drops previous pairs"
    nominal.sort(key=lambda r: (screen_order[r["screen_id"]], float(r["b_p"]), r["gene_id"]))

    context = {"vasc": "vascular", "neun": "neuronal", "neuropil": "neuropil"}
    columns = {
        "screen": "screen_id", "mediator": "mediator_symbol",
        "outcome": "gene_name", "outcome_gene_id": "gene_id",
        "P.Value_DE": "historical_p", "fdr_DE": "historical_q",
        "logFC_base": "c_beta", "P.Value_base": "c_p",
        "logFC_med": "cprime_beta", "P.Value_med": "cprime_p",
        "logFC_med_vec": "b_beta", "P.Value_med_vec": "b_p",
        "fdr_med_vec": "b_q", "fdr_med_vec_pooled": "b_q_global",
        "mediator_logFC_matched": "a_beta", "mediator_P.Value_matched": "a_p",
        "n_observations": "n_samples", "n_donors": "n_donors",
        "coefficient_shrinkage": "coefficient_shrinkage", "direction_reversal": "direction_reversal",
    }
    table = []
    for r in nominal:
        out = {name: r[field] for name, field in columns.items()}
        out["mediator_microenv"] = context[r["source"]]
        out["outcome_microenv"] = context[r["target"]]
        out["shrinkage_pct"] = 100 * (1 - float(r["cprime_beta"]) / float(r["c_beta"]))
        out["passed_previous_FDR10_screen"] = "TRUE" if key(r) in previous_keys else "FALSE"
        table.append(out)
    fields = ["screen", "mediator", "mediator_microenv", "outcome_microenv", "outcome", "outcome_gene_id",
              *list(columns)[4:10], "shrinkage_pct", *list(columns)[10:], "passed_previous_FDR10_screen"]
    counts_old = Counter(r["screen_id"] for r in previous)
    counts_new = Counter(r["screen_id"] for r in nominal)
    summary = [{"screen": sid, "previous_FDR10_pairs": counts_old[sid],
                "nominal_p05_pairs": counts_new[sid], "additional_pairs": counts_new[sid] - counts_old[sid],
                "unique_outcome_genes": len({r["gene_id"] for r in nominal if r["screen_id"] == sid})}
               for sid in screen_order]
    outdir = run / "report"
    outdir.mkdir(exist_ok=True)
    for extension, delimiter in (("csv", ","), ("tab", "\t")):
        write_table(outdir / f"mediation_hits_nominal_p05.{extension}", table, fields, delimiter)
    write_table(outdir / "nominal_screen_summary.tsv", summary, list(summary[0]))
    provenance = {
        "source": "primary/all_pairs.tsv.gz", "source_sha256": source_sha,
        "source_run_signature": next(iter(signatures)), "tested_pairs": len(rows),
        "filter": {"historical_p": "< 0.05", "c_p": "< 0.05", "b_p": "< 0.05", "cprime_p": ">= 0.05"},
        "matched_mediator_diagnosis_gate": False, "FDR_columns": "unchanged annotations; original full testing families",
        "previous_pairs": len(previous), "nominal_pairs": len(nominal),
        "unique_outcome_genes": len({r["gene_id"] for r in nominal}),
        "coefficient_shrinkage_pairs": sum(r["coefficient_shrinkage"] == "TRUE" for r in nominal),
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    }
    (outdir / "nominal_filter_provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    assert hashlib.sha256(source.read_bytes()).hexdigest() == source_sha, "Source results changed"
    print(json.dumps({"summary": summary, "previous_pairs": len(previous), "nominal_pairs": len(nominal),
                      "unique_outcome_genes": provenance["unique_outcome_genes"]}, indent=2))


if __name__ == "__main__":
    main()
