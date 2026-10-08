#!/usr/bin/env python3
## independently verify exported identifiers, model comparisons, and BH families.
import argparse
import csv
import gzip
import math
from pathlib import Path


def read_table(path):
    opener = gzip.open if str(path).endswith('.gz') else open
    with opener(path, 'rt', newline='') as handle:
        return list(csv.DictReader(handle, delimiter='\t'))


def bh(pvalues):
    n = len(pvalues)
    result = [None] * n
    previous = 1.0
    for rank, index in reversed(list(enumerate(sorted(range(n), key=pvalues.__getitem__), 1))):
        previous = min(previous, pvalues[index] * n / rank)
        result[index] = previous
    return result


def check_family(rows, prefix):
    pvalues = [float(row[prefix + '_p']) for row in rows]
    observed = [float(row[prefix + '_q']) for row in rows]
    assert all(math.isfinite(p) and 0 <= p <= 1 for p in pvalues)
    assert all(math.isclose(a, b, rel_tol=1e-10, abs_tol=1e-12)
               for a, b in zip(observed, bh(pvalues))), prefix + ' BH family mismatch'


def check_robustness(outdir, manifest, primary, signature):
    expected = set()
    status = read_table(outdir / 'sensitivity/status.tsv')
    voom_engine = any(row['analysis'] == 'fixed_weights' for row in status)
    for screen in manifest:
        labels = ['slide_rin'] + (['matched_logcounts', 'fixed_weights'] if voom_engine else ['voom'])
        if 'neun' in (screen['source'], screen['target']):
            labels.append('without_spd07')
        expected.update((screen['screen_id'], label) for label in labels)
    completed = {(row['screen_id'], row['analysis']) for row in status
                 if row['status'] == 'complete'}
    assert expected <= completed, 'Required sensitivity tasks missing or failed'
    assert all(row['status'] == 'complete' for row in status), 'Sensitivity task failed'
    assert {row['run_signature'] for row in status} == {signature}, 'Stale sensitivity status'
    for sid, label in sorted(expected):
        rows = read_table(outdir / 'sensitivity' / (sid + '_' + label + '.tsv.gz'))
        assert rows, 'Empty sensitivity table: ' + sid + '/' + label
        assert len({r['gene_id'] for r in rows}) == len(rows), 'Duplicate sensitivity gene'
        genes = {r['gene_id'] for r in rows}
        primary_genes = {r['gene_id'] for r in primary[sid]}
        if label in ('without_spd07', 'voom'):
            assert genes <= primary_genes, 'Refiltered universe exceeds primary universe'
        else:
            assert genes == primary_genes, \
                'Sensitivity testing universe differs from primary: ' + sid + '/' + label
        assert {r['run_signature'] for r in rows} == {signature}, 'Stale sensitivity table'
        assert {r['screen_id'] for r in rows} == {sid}
        assert {r['analysis'] for r in rows} == {label}
        for prefix in ('c', 'cprime', 'b'):
            check_family(rows, prefix)
        if label != 'without_spd07':
            for field in ('n_samples', 'n_donors'):
                assert {r[field] for r in rows} == {r[field] for r in primary[sid]}
        else:
            assert 0 < int(rows[0]['n_samples']) < int(primary[sid][0]['n_samples'])
            assert 0 < int(rows[0]['n_donors']) <= int(primary[sid][0]['n_donors'])
    overlap = read_table(outdir / 'overlap/status.tsv')
    assert {r['screen_id'] for r in overlap} == {s['screen_id'] for s in manifest}
    assert {r['run_signature'] for r in overlap} == {signature}, 'Stale overlap status'
    allowed = {'complete', 'not_required_no_primary_hits', 'untestable_or_failed'}
    assert all(r['status'] in allowed for r in overlap), 'Unknown overlap status'
    limitations = sum(r['status'] == 'untestable_or_failed' for r in overlap)
    print('Validated {} required sensitivity tables: complete tasks, current signatures, '
          'gene universes, sample counts, and full-gene BH corrections.'.format(len(expected)))
    print('Overlap status verified; {} explicitly untestable or failed refits.'.format(limitations))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--manifest', type=Path, default=Path(__file__).resolve().parents[1] / 'screens.tsv')
    parser.add_argument('--allow-partial', action='store_true')
    parser.add_argument('--require-robustness', action='store_true')
    args = parser.parse_args()
    assert not (args.allow_partial and args.require_robustness), 'Robustness requires a complete run'
    manifest = read_table(args.manifest)
    audit = {r['screen_id']: r for r in read_table(args.outdir / 'audit/matched_designs.tsv')}
    observed = {}
    for screen in manifest:
        sid = screen['screen_id']
        path = args.outdir / 'primary' / (sid + '_results.tsv.gz')
        if not path.exists():
            assert args.allow_partial, 'Missing screen: ' + sid
            continue
        rows = read_table(path)
        assert len(rows) == int(audit[sid]['target_genes'])
        assert len({r['gene_id'] for r in rows}) == len(rows)
        assert {r['screen_id'] for r in rows} == {sid}
        assert {int(r['n_samples']) for r in rows} == {int(screen['expected_samples'])}
        assert {int(r['n_donors']) for r in rows} == {int(screen['expected_donors'])}
        for prefix in ('c', 'cprime', 'b'):
            check_family(rows, prefix)
        samples = read_table(args.outdir / 'primary' / (sid + '_samples.tsv'))
        assert len(samples) == int(screen['expected_samples'])
        assert len({r['key'] for r in samples}) == len(samples)
        assert len({r['donor'] for r in samples}) == int(screen['expected_donors'])
        observed[sid] = rows
    assert observed, 'No completed screen tables'
    a = 'fgf1_neuropil_vasc'
    b = 'fgf2_neuropil_vasc'
    if a in observed and b in observed:
        first = {r['gene_id']: r for r in observed[a]}
        second = {r['gene_id']: r for r in observed[b]}
        assert first.keys() == second.keys()
        for gene in first:
            for field in ('c_beta', 'c_se', 'c_p', 'c_q'):
                assert first[gene][field] == second[gene][field], 'Shared baseline mismatch'
    combined = args.outdir / 'primary/all_pairs.tsv.gz'
    if len(observed) == len(manifest):
        assert combined.exists(), 'Combined testing family unavailable'
        rows = read_table(combined)
        assert len(rows) == sum(len(r) for r in observed.values())
        assert len({(r['screen_id'], r['gene_id']) for r in rows}) == len(rows)
        assert len({r['run_signature'] for r in rows}) == 1
        assert all(r['pooled_family_complete'] == 'TRUE' for r in rows)
        expected = bh([float(r['b_p']) for r in rows])
        assert all(math.isclose(float(r['b_q_global']), q, rel_tol=1e-10, abs_tol=1e-12)
                   for r, q in zip(rows, expected)), 'Pooled BH family mismatch'
    print('Validated {} screens / {} pairs: unique sample and gene IDs, audit counts, '
          'full-gene BH corrections, and identical shared baselines.'.format(
              len(observed), sum(len(r) for r in observed.values())))
    if len(observed) == len(manifest):
        print('Pooled BH correction across all five screens also verified.')
        if args.require_robustness:
            check_robustness(args.outdir, manifest, observed, rows[0]['run_signature'])


if __name__ == '__main__':
    main()
