"""Verify candidate plotting coverage, inputs, PDF pages and saved provenance."""
import argparse
import csv
import gzip
import hashlib
import json
import math
import subprocess
from collections import defaultdict
from pathlib import Path


def read(path, delimiter='\t'):
    opener = gzip.open if path.suffix == '.gz' else open
    with opener(path, 'rt') as f:
        return list(csv.DictReader(f, delimiter=delimiter))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def verify_pdf(path, expected):
    info = subprocess.check_output(['pdfinfo', str(path)], text=True)
    pages = int(next(line.split(':', 1)[1] for line in info.splitlines() if line.startswith('Pages:')))
    assert pages == len(expected), (path, pages)
    content = subprocess.check_output(['pdftotext', '-layout', str(path), '-'], text=True).split('\f')
    assert len(content) == pages + 1 and not content[-1].strip()
    for page, record in zip(content, expected):
        assert record['outcome_gene_id'] in page and record['screen'] in page
        assert record['outcome'] in page and record['mediator'] in page
        assert 'Outcome: baseline' in page and 'Outcome: + mediator' in page and 'Mediator: source context' in page
        assert 'BH FDR' in page and 'no new statistical tests' in page


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--source-run', required=True, type=Path)
    p.add_argument('--plots', required=True, type=Path)
    args = p.parse_args()
    run, plots = args.source_run, args.plots
    out = run / 'candidate_boxplots'
    hits = read(run / 'report/mediation_hits_nominal_p05.csv', ',')
    index = read(out / 'pair_plot_index.tsv')
    selected = read(out / 'selected_pairs.tsv')
    assert len(index) == len(hits) == 110
    assert len(selected) == 16 and len({r['outcome'] for r in selected}) == 11
    for page, (entry, hit) in enumerate(zip(index, hits), 1):
        assert int(entry['page_all']) == page
        assert all(entry[k] == v for k, v in hit.items())
    grouped = defaultdict(list)
    for row in read(out / 'plot_values.tsv.gz'):
        grouped[row['screen'], row['gene_id']].append(row)
        assert row['donor'].startswith('D') and len(row['donor']) == 4
        assert row['Dx'] in ('SCZ', 'NTC')
        for field in ('outcome_baseline', 'outcome_adjusted', 'mediator_adjusted'):
            assert math.isfinite(float(row[field]))
    assert set(grouped) == {(r['screen'], r['outcome_gene_id']) for r in hits}
    source_values = {}
    for hit in hits:
        rows = grouped[hit['screen'], hit['outcome_gene_id']]
        assert len(rows) == int(hit['n_observations'])
        assert len({r['donor'] for r in rows}) == int(hit['n_donors'])
        assert len({(r['donor'], r['SpD']) for r in rows}) == len(rows)
        assert abs(sum(float(r['outcome_baseline']) - float(r['outcome_adjusted']) for r in rows) / len(rows)) < 1e-12
        for r in rows:
            key = hit['screen'], r['donor'], r['SpD']
            value = r['Dx'], r['mediator_adjusted']
            assert source_values.setdefault(key, value) == value
    for filename in ('candidate_boxplots_all_110.pdf', 'candidate_boxplots_all_110_donor_means.pdf'):
        verify_pdf(plots / filename, hits)
    for filename in ('candidate_boxplots_selected.pdf', 'candidate_boxplots_selected_donor_means.pdf'):
        verify_pdf(plots / filename, selected)
    for row in selected:
        verify_pdf(plots / row['plot_file'], [row])
    provenance = json.loads((out / 'provenance.json').read_text())
    for filename, sha in provenance['inputs_sha256'].items():
        assert digest(run / filename) == sha
    for filename, sha in provenance['outputs_sha256'].items():
        assert digest(plots / filename) == sha
    for name, key in [('plot_candidates.py', 'script_sha256'), ('export_candidate_boxplots.R', 'export_script_sha256'), ('boxplot_shortlist.tsv', 'shortlist_sha256')]:
        assert digest(Path(__file__).resolve().parents[1] / name) == provenance[key]
    checks = read(out / 'model_verification.tsv')
    assert len(checks) == 5 and max(float(r['max_GLS_coefficient_error']) for r in checks) < 1e-8
    print('Verified 110 pairs, 16 selected pairs, 252 atlas pages plus 16 individual PDFs, plotting coverage, statistics index and file hashes.')


if __name__ == '__main__':
    main()
