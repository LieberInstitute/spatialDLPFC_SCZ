"""Generate gene-overlap and GO figures for the fixed nominal mediation screen."""

import argparse
import csv
import gzip
import hashlib
import json
import sqlite3
import shutil
import io
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import scipy
from scipy.stats import hypergeom


LABELS = {
    "ptn_vasc_neuropil": "PTN vascular → neuropil",
    "ptn_vasc_neun": "PTN vascular → neuronal",
    "fgf1_neuropil_vasc": "FGF1 neuropil → vascular",
    "fgf2_neuropil_vasc": "FGF2 neuropil → vascular",
    "fgf1_neun_vasc": "FGF1 neuronal → vascular",
}
ONTOLOGIES = {"CC": "Cellular component", "BP": "Biological process", "MF": "Molecular function"}


def read_table(path, delimiter="\t"):
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def write_table(path, rows, fields):
    buffer = io.StringIO(newline="")
    writer = csv.DictWriter(buffer, fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    data = buffer.getvalue().encode()
    path.write_bytes(gzip.compress(data, mtime=0) if path.suffix == ".gz" else data)


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def bh(values):
    values = np.asarray(values)
    order = np.argsort(values)
    result = np.empty(len(values))
    result[order] = np.minimum.accumulate((values[order] * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    return np.minimum(result, 1)


def save_figure(fig, basename):
    for extension in ("pdf", "svg", "png"):
        fig.savefig(basename.with_suffix("." + extension), dpi=300, facecolor="white")
    plt.close(fig)


def overlap_figure(hits, screens, out, plots, draw=True):
    sets = {s: {r["outcome_gene_id"] for r in hits if r["screen"] == s} for s in screens}
    names = {r["outcome_gene_id"]: r["outcome"] for r in hits}
    genes = sorted(set.union(*sets.values()))
    membership = [{"gene_id": g, "gene_symbol": names[g], **{s: int(g in sets[s]) for s in screens}} for g in genes]
    write_table(out / "gene_overlap_membership.tsv", membership, ["gene_id", "gene_symbol", *screens])
    groups = defaultdict(list)
    for g in genes:
        groups[tuple(s for s in screens if g in sets[s])].append(g)
    groups = sorted(groups.items(), key=lambda x: (-len(x[1]), x[0]))
    intersections = [{"screens": ";".join(k), "n_genes": len(v), "gene_ids": ";".join(v),
                      "gene_symbols": ";".join(names[g] for g in v)} for k, v in groups]
    write_table(out / "gene_overlap_intersections.tsv", intersections, list(intersections[0]))
    assert sum(len(v) for _, v in groups) == len(genes)
    assert sum(len(k) * len(v) for k, v in groups) == len(hits)

    if not draw:
        return
    fig = plt.figure(figsize=(12, 6.1))
    grid = fig.add_gridspec(2, 2, width_ratios=[1.45, 3.8], height_ratios=[2.4, 1.6],
                           left=.245, right=.97, bottom=.16, top=.81, hspace=.06, wspace=.07)
    top = fig.add_subplot(grid[0, 1])
    matrix = fig.add_subplot(grid[1, 1], sharex=top)
    bars = fig.add_subplot(grid[1, 0], sharey=matrix)
    x = np.arange(len(groups))
    counts = [len(v) for _, v in groups]
    colors = ["#24788A" if len(k) > 1 else "#737B83" for k, _ in groups]
    top.bar(x, counts, color=colors, width=.72)
    for i, n in enumerate(counts):
        top.text(i, n + .6, str(n), ha="center", va="bottom", fontsize=11)
    top.set_ylim(0, max(counts) * 1.16)
    top.set_ylabel("Genes in exact intersection")
    top.tick_params(axis="x", bottom=False, labelbottom=False)
    top.spines[["top", "right", "bottom"]].set_visible(False)
    top.grid(axis="y", color="#E8EBEE", linewidth=.6)
    top.set_axisbelow(True)
    for i, (members, _) in enumerate(groups):
        matrix.scatter([i] * len(screens), range(len(screens)), s=38, color="#E0E4E7", zorder=1)
        ys = [j for j, s in enumerate(screens) if s in members]
        if len(ys) > 1:
            matrix.plot([i, i], [min(ys), max(ys)], color=colors[i], linewidth=2, zorder=2)
        matrix.scatter([i] * len(ys), ys, s=48, color=colors[i], zorder=3)
    matrix.set_ylim(len(screens) - .6, -.6)
    matrix.set_xlim(-.65, len(groups) - .35)
    matrix.set_xticks([])
    matrix.tick_params(axis="y", left=False, labelleft=False)
    matrix.spines[:].set_visible(False)
    bars.barh(range(len(screens)), [len(sets[s]) for s in screens], height=.6, color="#24788A")
    bars.set_yticks(range(len(screens)), [LABELS[s] for s in screens], fontsize=10)
    bars.set_xlim(max(map(len, sets.values())) * 1.2, 0)
    bars.set_xlabel("Set size")
    bars.spines[["top", "left", "right"]].set_visible(False)
    bars.tick_params(axis="y", length=0)
    for i, s in enumerate(screens):
        bars.text(len(sets[s]) + 1, i, str(len(sets[s])), ha="right", va="center", fontsize=10)
    fig.text(.03, .94, "B", fontsize=23, weight="bold")
    fig.text(.085, .945, "Outcome-gene overlap across nominal mediation screens", fontsize=15, weight="bold")
    fig.text(.085, .897, f"{len(hits)} mediator–outcome pairs · {len(genes)} unique outcome genes · 5 screens", color="#4D5964", fontsize=11)
    fig.text(.245, .045, "Intersections use gene IDs. The same gene in different target microenvironments counts once.\n"
             "Teal marks overlap between screens; selection uses nominal mediation thresholds.", fontsize=9, color="#4D5964")
    save_figure(fig, plots / "figure_6B_gene_overlap")


def annotations(org_path, go_path, gene_ids, out):
    org = sqlite3.connect(f"file:{org_path}?mode=ro", uri=True)
    go = sqlite3.connect(f"file:{go_path}?mode=ro", uri=True)
    mapping = defaultdict(set)
    for ensembl, entrez in org.execute("SELECT DISTINCT e.ensembl_id,g.gene_id FROM ensembl e JOIN genes g USING(_id)"):
        if ensembl in gene_ids:
            mapping[ensembl].add(entrez)
    symbols = dict(org.execute("SELECT g.gene_id,i.symbol FROM genes g JOIN gene_info i USING(_id)"))
    mapped_genes = set.union(set(), *mapping.values())
    org.execute("CREATE TEMP TABLE selected (entrez TEXT PRIMARY KEY)")
    org.executemany("INSERT INTO selected VALUES (?)", [(g,) for g in sorted(mapped_genes)])
    terms = {term: (name, ont) for term, name, ont in go.execute("SELECT go_id,term,ontology FROM go_term")}
    term_genes = {}
    snapshots = []
    for ont in ONTOLOGIES:
        sets = defaultdict(set)
        query = f"SELECT DISTINCT a.go_id,g.gene_id FROM go_{ont.lower()}_all a JOIN genes g USING(_id) JOIN selected s ON s.entrez=g.gene_id"
        for term, gene in org.execute(query):
            if term in terms and terms[term][1] == ont:
                sets[term].add(gene)
                snapshots.append({"ontology": ont, "go_id": term, "entrez_id": gene})
        term_genes[ont] = sets
    snapshot_fields = ["ontology", "go_id", "entrez_id"]
    snapshots.sort(key=lambda r: tuple(r[f] for f in snapshot_fields))
    write_table(out / "GO_annotation_snapshot.tsv.gz", snapshots, snapshot_fields)
    write_table(out / "GO_gene_mapping.tsv", [{"gene_id": g, "entrez_id": e, "symbol": symbols[e]}
                for g in sorted(mapping) for e in sorted(mapping[g])], ["gene_id", "entrez_id", "symbol"])
    used_terms = set.union(*(set(x) for x in term_genes.values()))
    write_table(out / "GO_term_snapshot.tsv", [{"go_id": t, "description": terms[t][0], "ontology": terms[t][1]}
                for t in sorted(used_terms)], ["go_id", "description", "ontology"])
    metadata = {"org.Hs.eg.db": dict(org.execute("SELECT * FROM metadata")), "GO.db": dict(go.execute("SELECT * FROM metadata"))}
    org.close()
    go.close()
    return mapping, symbols, terms, term_genes, metadata


def saved_annotations(directory):
    mapping, symbols = defaultdict(set), {}
    for r in read_table(directory / "GO_gene_mapping.tsv"):
        mapping[r["gene_id"]].add(r["entrez_id"])
        symbols[r["entrez_id"]] = r["symbol"]
    terms = {r["go_id"]: (r["description"], r["ontology"]) for r in read_table(directory / "GO_term_snapshot.tsv")}
    term_genes = {ont: defaultdict(set) for ont in ONTOLOGIES}
    for r in read_table(directory / "GO_annotation_snapshot.tsv.gz"):
        term_genes[r["ontology"]][r["go_id"]].add(r["entrez_id"])
    provenance = json.loads((directory / "provenance.json").read_text())
    return mapping, symbols, terms, term_genes, provenance["annotation_metadata"]


def background_genes(rows, screen, mode):
    return {r["gene_id"] for r in rows if r["screen_id"] == screen
            and (mode == "all_tested" or (float(r["historical_p"]) < .05 and float(r["c_p"]) < .05))}


def background_caption(mode):
    return ("Background: all tested outcome genes in each target microenvironment."
            if mode == "all_tested" else
            "Background: outcome genes passing manuscript and matched-baseline diagnosis p < 0.05 in that screen.")


def go_results(rows, hits, screens, mapping, symbols, terms, term_genes, out, mode):
    results, summaries, universe_rows = [], [], []
    for screen in screens:
        eligible = background_genes(rows, screen, mode)
        selected = {r["outcome_gene_id"] for r in hits if r["screen"] == screen}
        assert selected <= eligible
        universe = set.union(set(), *(mapping[g] for g in eligible))
        candidates = set.union(set(), *(mapping[g] for g in selected))
        universe_rows.extend({"screen": screen, "gene_id": g, "selected": int(g in selected),
                              "entrez_ids": ";".join(sorted(mapping[g]))} for g in sorted(eligible))
        for ont in ONTOLOGIES:
            annotated = set.union(set(), *term_genes[ont].values())
            background = universe & annotated
            query_genes = candidates & background
            N, n = len(background), len(query_genes)
            family = []
            all_term_count = 0
            for term, annotated_genes in sorted(term_genes[ont].items()):
                members = annotated_genes & background
                overlap = members & query_genes
                ## match the ERC enrichGO/DOSE family: query-overlapping terms of eligible size.
                if not 10 <= len(members) <= 500:
                    continue
                all_term_count += 1
                if not overlap:
                    continue
                k, M = len(overlap), len(members)
                family.append({"screen": screen, "ontology": ont, "go_id": term, "description": terms[term][0],
                               "Count": k, "query_annotated": n, "term_background": M, "background_annotated": N,
                               "GeneRatio": k / n, "BgRatio": M / N,
                               "fold_enrichment": (k / n) / (M / N), "pvalue": float(hypergeom.sf(k - 1, N, M, n)),
                               "entrez_ids": ";".join(sorted(overlap)),
                               "gene_symbols": ";".join(symbols[g] for g in sorted(overlap))})
            for row, adjusted in zip(family, bh([r["pvalue"] for r in family])):
                row["p_adjust"] = float(adjusted)
            ## sensitivity: include zero-overlap terms (p = 1) in each BH family.
            expanded = bh([r["pvalue"] for r in family] + [1.0] * (all_term_count - len(family)))
            for row, adjusted in zip(family, expanded):
                row["p_adjust_all_eligible_terms"] = float(adjusted)
            results.extend(family)
            summaries.append({"screen": screen, "ontology": ont, "background_mode": mode, "background_ensembl": len(eligible),
                              "mapped_universe_entrez": len(universe), "background_annotated": N,
                              "candidate_ensembl": len(selected), "candidate_unmapped": sum(not mapping[g] for g in selected),
                              "query_annotated": n, "tested_terms": len(family), "all_eligible_terms": all_term_count,
                              "BH05_terms": sum(r["p_adjust"] < .05 for r in family),
                              "BH05_terms_2plus_genes": sum(r["p_adjust"] < .05 and r["Count"] >= 2 for r in family)})
    ## an additional correction across all screen/ontology tests is reported for context.
    for r, adjusted in zip(results, bh([r["pvalue"] for r in results])):
        r["p_adjust_all_screens_ontologies"] = float(adjusted)
    fields = list(results[0])
    results.sort(key=lambda r: (screens.index(r["screen"]), r["ontology"], r["pvalue"], r["go_id"]))
    write_table(out / "GO_results_all.tsv.gz", results, fields)
    write_table(out / "GO_significant_BH05.tsv", [r for r in results if r["p_adjust"] < .05], fields)
    write_table(out / "GO_group_summary.tsv", summaries, list(summaries[0]))
    write_table(out / "GO_eligible_universes.tsv", universe_rows, list(universe_rows[0]))
    return results, summaries


def go_figure(results, screens, ont, plots, mode):
    significant = [r for r in results if r["ontology"] == ont and r["p_adjust"] < .05 and r["Count"] >= 2]
    display_terms = set()
    for screen in screens:
        group = sorted((r for r in significant if r["screen"] == screen), key=lambda r: (r["p_adjust"], r["pvalue"], r["go_id"]))
        display_terms.update(r["go_id"] for r in group[:3])
    ordered = sorted(display_terms, key=lambda term: (min(r["p_adjust"] for r in significant if r["go_id"] == term), term))
    shown = [r for r in significant if r["go_id"] in display_terms]
    height = max(5.6, len(ordered) * .65 + 2.2)
    fig, ax = plt.subplots(figsize=(12, height))
    fig.subplots_adjust(left=.39, right=.84, top=.84, bottom=.27)
    if shown:
        xs = [screens.index(r["screen"]) for r in shown]
        ys = [ordered.index(r["go_id"]) for r in shown]
        points = ax.scatter(xs, ys, s=[40 + 650 * r["GeneRatio"] for r in shown],
                            c=[-np.log10(max(r["p_adjust"], 1e-300)) for r in shown],
                            cmap="viridis_r", vmin=-np.log10(.05), edgecolors="#3D4852", linewidths=.5)
        ax.set_yticks(range(len(ordered)), [textwrap.fill(next(r["description"] for r in shown if r["go_id"] == t), 43) for t in ordered])
        ax.set_ylim(len(ordered) - .6, -.6)
        colorbar = fig.colorbar(points, ax=ax, fraction=.045, pad=.065)
        colorbar.set_label("−log10 GO BH-adjusted p", fontsize=10)
        for ratio in (.1, .25, .5):
            ax.scatter([], [], s=40 + 650 * ratio, label=f"{ratio:.2f}", color="#C5CCD2", edgecolors="#3D4852", linewidths=.5)
        ax.legend(title="Gene ratio", loc="upper left", bbox_to_anchor=(1.26, .55), frameon=False, labelspacing=1.25)
    else:
        ax.text(.5, .5, "No GO terms meet BH-adjusted p < 0.05\nwith at least two candidate genes.", ha="center", va="center", transform=ax.transAxes, fontsize=12)
        ax.set_yticks([])
    ax.set_xticks(range(len(screens)), [LABELS[s].replace(" → ", "\n→ ") for s in screens], rotation=45, ha="right", fontsize=9)
    ax.set_xlim(-.6, len(screens) - .4)
    ax.grid(color="#E7EBEE", linewidth=.6)
    ax.set_axisbelow(True)
    ax.spines[["top", "right"]].set_visible(False)
    fig.text(.03, .94, "C", fontsize=23, weight="bold")
    fig.text(.085, .945, f"GO enrichment · {ONTOLOGIES[ont].lower()}", fontsize=15, weight="bold")
    fig.text(.085, .895, "Nominal mediation candidates; GO p-values use BH correction within each screen and ontology.", fontsize=10, color="#4D5964")
    fig.text(.085, .035, background_caption(mode) + "\n"
             "Top three qualifying terms per screen; displayed dots require GO BH p < 0.05 and ≥2 candidate genes.\n"
             "Separate backgrounds and annotation coverage can change gene-ratio denominators between screens.", fontsize=9, color="#4D5964")
    save_figure(fig, plots / f"figure_6C_GO_{ont}")


def go_summary_figure(results, screens, plots, mode):
    onts = list(ONTOLOGIES)
    minimum = np.array([[min((r["p_adjust"] for r in results if r["screen"] == s and r["ontology"] == ont), default=1)
                         for ont in onts] for s in screens])
    fig, ax = plt.subplots(figsize=(10, 5.7))
    fig.subplots_adjust(left=.34, right=.88, top=.77, bottom=.24)
    image = ax.imshow(-np.log10(minimum), cmap="Blues", vmin=0, vmax=max(-np.log10(.05), float(-np.log10(minimum.min()))), aspect="auto")
    for i in range(len(screens)):
        for j in range(len(onts)):
            ax.text(j, i, f"{minimum[i, j]:.3g}", ha="center", va="center",
                    color="white" if -np.log10(minimum[i, j]) > .5 * max(-np.log10(.05), float(-np.log10(minimum.min()))) else "#253746", fontsize=12)
    ax.set_yticks(range(len(screens)), [LABELS[s] for s in screens], fontsize=10)
    ax.set_xticks(range(len(onts)), [ONTOLOGIES[o].replace(" ", "\n") for o in onts], fontsize=10)
    ax.tick_params(length=0, pad=9)
    ax.spines[:].set_visible(False)
    colorbar = fig.colorbar(image, ax=ax, fraction=.045, pad=.05)
    colorbar.set_label("−log10 minimum BH p", fontsize=10)
    fig.text(.03, .94, "C", fontsize=23, weight="bold")
    fig.text(.09, .945, "GO enrichment across nominal mediation screens", fontsize=14, weight="bold")
    n_significant = sum(r["p_adjust"] < .05 for r in results)
    message = ("No terms pass BH-adjusted p < 0.05 in any screen or ontology." if not n_significant
               else f"{n_significant} screen/ontology–term associations pass BH-adjusted p < 0.05.")
    fig.text(.09, .88, message, fontsize=11, color="#253746")
    fig.text(.09, .835, "Cells show the smallest GO BH-adjusted p in each testing family.", fontsize=10, color="#4D5964")
    fig.text(.09, .065, background_caption(mode) + "\n"
             "Hypergeometric tests; inherited GO annotations; eligible term sizes 10–500 background genes.\n"
             f"BH correction within each screen/ontology; {len(results):,} term tests in total. Values are adjusted p-values.", fontsize=9, color="#4D5964")
    save_figure(fig, plots / "figure_6C_GO_summary")


def background_comparison(run, plots):
    directories = [run / "downstream_nominal", run / "downstream_nominal_all_tested"]
    provenance = [json.loads((d / "provenance.json").read_text()) for d in directories]
    assert provenance[0]["input_sha256"] == provenance[1]["input_sha256"]
    matrices, counts, rows = [], [], []
    for d, mode in zip(directories, ("eligible", "all_tested")):
        results = read_table(d / "GO_results_all.tsv.gz")
        matrix = []
        for screen in LABELS:
            line = []
            for ont in ONTOLOGIES:
                family = [r for r in results if r["screen"] == screen and r["ontology"] == ont]
                minimum = min(float(r["p_adjust"]) for r in family)
                line.append(minimum)
                rows.append({"background": mode, "screen": screen, "ontology": ont,
                             "minimum_BH_p": minimum,
                             "BH05_terms": sum(float(r["p_adjust"]) < .05 for r in family),
                             "pooled_BH05_terms": sum(float(r["p_adjust_all_screens_ontologies"]) < .05 for r in family),
                             "expanded_family_BH05_terms": sum(float(r["p_adjust_all_eligible_terms"]) < .05 for r in family)})
            matrix.append(line)
        matrices.append(np.array(matrix))
        counts.append(sum(float(r["p_adjust"]) < .05 for r in results))
    write_table(run / "report/GO_background_comparison.tsv", rows, list(rows[0]))
    fig, axes = plt.subplots(1, 2, figsize=(12, 5.6), sharey=True)
    fig.subplots_adjust(left=.28, right=.89, top=.73, bottom=.23, wspace=.12)
    vmax = max(float(-np.log10(m.min())) for m in matrices)
    for ax, m, title, count in zip(axes, matrices,
                                  ("Screen-eligible DEGs", "All tested target genes"), counts):
        im = ax.imshow(-np.log10(m), cmap="Blues", vmin=0, vmax=vmax, aspect="auto")
        for i in range(len(LABELS)):
            for j in range(len(ONTOLOGIES)):
                ax.text(j, i, f"{m[i, j]:.3g}", ha="center", va="center", fontsize=11,
                        color="white" if -np.log10(m[i, j]) > vmax / 2 else "#253746")
        ax.set_title(f"{title}\n{count} term associations at BH p < 0.05", fontsize=11, pad=13)
        ax.set_xticks(range(3), list(ONTOLOGIES))
        ax.set_yticks(range(5), list(LABELS.values()), fontsize=10)
        ax.tick_params(length=0, pad=9)
        ax.spines[:].set_visible(False)
    cax = fig.add_axes([.915, .23, .018, .5])
    fig.colorbar(im, cax=cax).set_label("−log10 minimum BH p", fontsize=10)
    fig.text(.03, .94, "C", fontsize=23, weight="bold")
    fig.text(.085, .945, "GO enrichment depends on the reference gene set", fontsize=15, weight="bold")
    fig.text(.085, .88, "Same 110 candidate pairs and GO annotations; cells show the smallest BH-adjusted p per family.", fontsize=10)
    fig.text(.085, .07, "CC: cellular component; BP: biological process; MF: molecular function. BH correction within each screen/ontology.\n"
             "Left: functional enrichment among eligible DEGs. Right: characterization relative to the tested target transcriptome.\n"
             "GO enrichment does not validate individual mediation pairs; overlapping terms can reflect the same genes.", fontsize=9, color="#4D5964")
    save_figure(fig, plots / "figure_6C_background_comparison")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-run", type=Path, required=True)
    parser.add_argument("--orgdb", type=Path)
    parser.add_argument("--godb", type=Path)
    parser.add_argument("--annotation-snapshot", type=Path)
    parser.add_argument("--plots", type=Path, required=True)
    parser.add_argument("--background", choices=("eligible", "all_tested"), default="eligible")
    parser.add_argument("--compare-backgrounds", action="store_true", help="Write comparison beside the eligible-background plots after both analyses exist")
    args = parser.parse_args()
    run, plots = args.source_run.resolve(), args.plots.resolve()
    out = run / ("downstream_nominal" if args.background == "eligible" else "downstream_nominal_all_tested")
    out.mkdir(exist_ok=True)
    plots.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "pdf.fonttype": 42, "svg.fonttype": "none"})
    source = run / "primary/all_pairs.tsv.gz"
    hit_file = run / "report/mediation_hits_nominal_p05.csv"
    fingerprints = {"source_results": sha256(source), "nominal_hits": sha256(hit_file)}
    if args.annotation_snapshot:
        snapshot_provenance = json.loads((args.annotation_snapshot / "provenance.json").read_text())
        assert snapshot_provenance["background_mode"] == args.background, "Snapshot background mismatch"
        for key in ("source_results", "nominal_hits"):
            assert snapshot_provenance["input_sha256"][key] == fingerprints[key], "Annotation snapshot belongs to different screening inputs"
        for filename in ("GO_gene_mapping.tsv", "GO_term_snapshot.tsv", "GO_annotation_snapshot.tsv.gz"):
            assert sha256(args.annotation_snapshot / filename) == snapshot_provenance["annotation_snapshot_sha256"][filename], "Annotation snapshot changed"
            if args.annotation_snapshot.resolve() != out:
                shutil.copy2(args.annotation_snapshot / filename, out / filename)
        fingerprints.update({k: v for k, v in snapshot_provenance["input_sha256"].items() if k.endswith(".sqlite")})
    else:
        assert args.orgdb and args.godb, "Supply both databases or --annotation-snapshot"
        fingerprints.update({"org.Hs.eg.sqlite": sha256(args.orgdb), "GO.sqlite": sha256(args.godb)})
    rows = read_table(source)
    hits = read_table(hit_file, ",")
    screens = list(LABELS)
    assert len({(r["screen_id"], r["gene_id"]) for r in rows}) == len(rows)
    assert {r["screen_id"] for r in rows} == set(screens)
    assert len({(r["screen"], r["outcome_gene_id"]) for r in hits}) == len(hits)
    expected = {(r["screen_id"], r["gene_id"]) for r in rows if float(r["historical_p"]) < .05
                and float(r["c_p"]) < .05 and float(r["b_p"]) < .05 and float(r["cprime_p"]) >= .05}
    assert expected == {(r["screen"], r["outcome_gene_id"]) for r in hits}
    overlap_figure(hits, screens, out, plots, draw=args.background == "eligible")
    eligible = set.union(*(background_genes(rows, screen, args.background) for screen in screens))
    if args.annotation_snapshot:
        mapping, symbols, terms, term_genes, metadata = saved_annotations(args.annotation_snapshot)
    else:
        mapping, symbols, terms, term_genes, metadata = annotations(args.orgdb, args.godb, eligible, out)
    results, summaries = go_results(rows, hits, screens, mapping, symbols, terms, term_genes, out, args.background)
    for ont in ONTOLOGIES:
        go_figure(results, screens, ont, plots, args.background)
    go_summary_figure(results, screens, plots, args.background)
    provenance = {"input_sha256": fingerprints, "annotation_packages": {"org.Hs.eg.db": "3.22.0", "GO.db": "3.22.0"},
                  "annotation_metadata": metadata, "software": {"numpy": np.__version__, "scipy": scipy.__version__, "matplotlib": matplotlib.__version__},
                  "script_sha256": sha256(Path(__file__)), "overlap_unit": "unique Ensembl outcome gene ID across all target contexts",
                  "background_mode": args.background, "GO_background": background_caption(args.background),
                  "GO_additional_BH": "all eligible terms including zero overlaps at p=1; separately pooled across the 15 query-overlap families",
                  "GO_qvalue": "Storey q-values are not used; significance uses BH p < .05",
                  "GO_test": "one-sided hypergeometric overrepresentation; inherited GO annotations; term background sizes 10–500",
                  "GO_BH_family": "all query-overlapping eligible terms within one screen and ontology, before plotting filters",
                  "GO_display": "top 3 BH<.05 terms with >=2 candidate genes per screen; union across screens",
                  "GO_sources": ["https://bioconductor.org/packages/3.22/data/annotation/html/org.Hs.eg.db.html", "https://bioconductor.org/packages/3.22/data/annotation/html/GO.db.html"],
                  "method_reference": "LFF_spatial_ERC/code/22_Mediation/06_Mediation_GO.R and DOSE 4.4.0 enricher_internal",
                  "mediation_pairs": len(hits), "unique_outcome_genes": len({r["outcome_gene_id"] for r in hits}),
                  "annotation_snapshot_sha256": {filename: sha256((args.annotation_snapshot or out) / filename)
                    for filename in ("GO_gene_mapping.tsv", "GO_term_snapshot.tsv", "GO_annotation_snapshot.tsv.gz")}}
    (out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
    if args.compare_backgrounds:
        background_comparison(run, plots.parent if args.background == "all_tested" else plots)
    print(json.dumps(summaries, indent=2))


if __name__ == "__main__":
    main()
