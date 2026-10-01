#!/usr/bin/env python3
## copy completed outputs into the analysis folder with portable provenance paths.
import argparse
import csv
import gzip
import hashlib
import re
from pathlib import Path


def portable_paths(text):
    text = re.sub(r'/(?:[^/\s]+/)*spatialDLPFC_SCZ(?:_LIBD4100)?(?=/|\s|$)',
                  '<PROJECT_ROOT>', text)
    return re.sub(r'/(?:Users|users|home)/[^/\s]+', '<HOME>', text)


def main():
    analysis = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-root', type=Path, default=analysis.parents[2])
    parser.add_argument('--label', required=True)
    args = parser.parse_args()
    assert re.fullmatch(r'[A-Za-z0-9_-]+', args.label), 'Invalid snapshot label'
    dest = analysis / 'reports' / args.label
    assert not dest.exists(), 'Snapshot already exists; use a new label'
    root = args.project_root.resolve()
    source = root / 'processed-data/17_mediation'
    plots = root / 'plots/17_mediation'
    assert (source / 'report/RESULTS.md').is_file(), 'Completed report unavailable'
    assert (source / 'sensitivity/status.tsv').is_file(), 'Sensitivity status unavailable'
    candidates = []
    for base, prefix in ((source, Path()), (plots, Path('figures'))):
        for path in sorted(base.rglob('*')):
            relative = path.relative_to(base)
            if not path.is_file() or path.name.startswith('.') or 'cache' in relative.parts:
                continue
            if path.suffix == '.rds' and 'provenance' not in relative.parts:
                continue
            candidates.append((path, prefix / relative))
    dest.mkdir(parents=True)
    manifest = []
    text_suffixes = {'.md', '.tsv', '.csv', '.txt', '.log'}
    for path, relative in candidates:
        original = path.read_bytes()
        archived = original
        normalized = False
        if path.suffix in text_suffixes or path.name.endswith(('.tsv.gz', '.csv.gz')):
            compressed = path.suffix == '.gz'
            payload = gzip.decompress(original) if compressed else original
            portable = portable_paths(payload.decode('utf-8')).encode('utf-8')
            normalized = portable != payload
            if normalized:
                archived = gzip.compress(portable, mtime=0) if compressed else portable
        target = dest / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        target.write_bytes(archived)
        manifest.append(dict(path=str(relative), source=str(path.relative_to(root)),
            source_sha256=hashlib.sha256(original).hexdigest(),
            archived_sha256=hashlib.sha256(archived).hexdigest(),
            source_bytes=len(original), archived_bytes=len(archived),
            path_normalized=normalized))
    with (dest / 'archive_manifest.tsv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=list(manifest[0]), delimiter='\t')
        writer.writeheader()
        writer.writerows(manifest)
    print('Archived {} files ({} path-normalized) in reports/{}.'.format(
        len(manifest), sum(r['path_normalized'] for r in manifest), args.label))


if __name__ == '__main__':
    main()
