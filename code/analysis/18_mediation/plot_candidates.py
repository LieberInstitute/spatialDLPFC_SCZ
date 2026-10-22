"""Plot all nominal candidate pairs and a documented biological shortlist."""

import argparse
import csv
import gzip
import hashlib
import json
import textwrap
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np

from downstream_nominal import read_table, sha256, write_table

CONTEXT = {'vasc': 'vascular', 'neun': 'neuronal', 'neuropil': 'neuropil'}
COLORS = {'NTC': '#3879A8', 'SCZ': '#C55C51'}
FIELDS = ('outcome_baseline', 'outcome_adjusted', 'mediator_adjusted')


def points_for_pair(rows, donor_means):
    if not donor_means:
        return rows
    by_donor = defaultdict(list)
    for row in rows:
        by_donor[row['donor']].append(row)
    result = []
    for donor, group in sorted(by_donor.items()):
        assert len({r['Dx'] for r in group}) == 1
        result.append({'donor': donor, 'Dx': group[0]['Dx'],
                       **{f: float(np.mean([float(r[f]) for r in group])) for f in FIELDS}})
    return result


def make_figure(hit, fitted, rows, donor_means=False, rationale=None):
    data = points_for_pair(rows, donor_means)
    fig, axes = plt.subplots(1, 3, figsize=(12.6, 6.2))
    fig.subplots_adjust(left=.075, right=.985, top=.69, bottom=.27, wspace=.34)
    target, source = CONTEXT[fitted['target']], CONTEXT[fitted['source']]
    fig.text(.055, .95, f"{hit['outcome']} [{target}]  |  {hit['mediator']} [{source}]", fontsize=18, weight='bold')
    unit = 'Donor means across available matched SpDs' if donor_means else 'Matched donor–SpD observations'
    fig.text(.055, .905, f"{unit}  ·  {hit['n_donors']} donors  ·  {hit['n_observations']} observations", fontsize=11, color='#455564')
    fig.text(.055, .867, f"{hit['screen']}  |  {hit['outcome_gene_id']}  |  Nominal screening candidate", fontsize=9, color='#586774')
    combined = [float(r[f]) for r in data for f in FIELDS[:2]]
    lim = (min(combined), max(combined))
    span = max(lim[1] - lim[0], .15)
    outcome_limits = (lim[0] - .1 * span, lim[1] + .12 * span)
    titles = ('Outcome: baseline', 'Outcome: + mediator', 'Mediator: source context')
    prefixes = ('c', 'cprime', 'a')
    seed = int(hashlib.sha256(hit['screen'].encode()).hexdigest()[:8], 16)
    rng = np.random.default_rng(seed)
    jitter = {r['donor'] + (r.get('SpD', '') if not donor_means else ''): rng.uniform(-.16, .16) for r in data}
    for panel, (ax, field, title, prefix) in enumerate(zip(axes, FIELDS, titles, prefixes)):
        ax.set_title(title, fontsize=12, weight='bold', pad=46,
                     bbox={'facecolor': '#DCEAF5' if panel == 2 else '#EAECEF', 'edgecolor': 'none', 'pad': 7})
        stats = f"ΔDx = {float(fitted[prefix + '_beta']):+.3f}   p = {float(fitted[prefix + '_p']):.3g}\nBH FDR = {float(fitted[prefix + '_q']):.3g}"
        ax.text(.5, 1.04, stats, ha='center', va='bottom', transform=ax.transAxes, fontsize=10)
        labels = []
        for x, dx in enumerate(COLORS):
            group = [r for r in data if r['Dx'] == dx]
            values = np.array([float(r[field]) for r in group])
            artist = ax.boxplot([values], positions=[x], widths=.52, patch_artist=True,
                showfliers=False, whis=1.5, manage_ticks=False,
                medianprops={'color': '#202E38', 'linewidth': 1.5},
                whiskerprops={'color': COLORS[dx]}, capprops={'color': COLORS[dx]})
            artist['boxes'][0].set(facecolor=COLORS[dx], edgecolor=COLORS[dx], alpha=.23)
            xs = [x + jitter[r['donor'] + (r.get('SpD', '') if not donor_means else '')] for r in group]
            ax.scatter(xs, values, s=13 if not donor_means else 23, color=COLORS[dx], alpha=.55, edgecolors='none', zorder=3)
            n_donors = len({r['donor'] for r in group})
            labels.append(f"{dx}\n{n_donors} donors" if donor_means else f"{dx}\n{len(group)} obs; {n_donors} donors")
        ax.set_xticks([0, 1], labels, fontsize=10)
        ax.set_xlim(-.5, 1.5)
        if panel < 2:
            ax.set_ylim(outcome_limits)
        ax.set_ylabel('Adjusted logcounts' if panel != 1 else '', fontsize=11)
        ax.spines[['top', 'right']].set_visible(False)
        ax.grid(axis='y', color='#E4E9EE', linewidth=.6)
        ax.set_axisbelow(True)
    fig.text(.055, .177, f"Conditional {hit['mediator']}–{hit['outcome']} association: b = {float(fitted['b_beta']):+.3f} per mediator SD, "
             f"p = {float(fitted['b_p']):.3g}, BH FDR = {float(fitted['b_q']):.3g}  |  Dx attenuation = {float(hit['shrinkage_pct']):.1f}%", fontsize=10)
    caption = ('Values remove centered age, sex and SpD effects; middle panel also removes the centered mediator effect.\n'
               'ΔDx = SCZ − NTC. P-values/FDR are from the original donor-blocked limma fits; boxes supply no new statistical tests.\n')
    caption += ('Each point averages adjusted values across the available SpDs of one donor; this is a descriptive display, not a refit.'
                if donor_means else 'Each point is a donor–SpD observation; multiple points from one donor are not independent. Shared outcome y-axis.')
    fig.text(.055, .082, caption, fontsize=8.7, linespacing=1.5, color='#455564')
    if rationale:
        fig.text(.055, .018, textwrap.fill('Selection: ' + rationale, width=165, max_lines=2, placeholder='…'), fontsize=8, color='#455564')
    else:
        fig.text(.055, .022, 'Exploratory association: significance loss after adjustment does not establish causal or complete mediation.', fontsize=9, color='#455564')
    fig.canvas.draw()
    bounds = fig.bbox
    for text in [*fig.texts, *(t for ax in axes for t in [ax.title, *ax.texts, *ax.get_xticklabels(), ax.yaxis.label])]:
        box = text.get_window_extent(fig.canvas.get_renderer())
        assert box.x0 >= bounds.x0 - 1 and box.y0 >= bounds.y0 - 1 and box.x1 <= bounds.x1 + 1 and box.y1 <= bounds.y1 + 1, text.get_text()
    return fig


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--source-run', type=Path, required=True)
    parser.add_argument('--plots', type=Path, required=True)
    args = parser.parse_args()
    run, plots = args.source_run, args.plots
    out = run / 'candidate_boxplots'
    plots.mkdir(parents=True, exist_ok=True)
    (plots / 'selected').mkdir(exist_ok=True)
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 10, 'pdf.fonttype': 42, 'svg.fonttype': 'none'})
    hitfile = run / 'report/mediation_hits_nominal_p05.csv'
    sourcefile = run / 'primary/all_pairs.tsv.gz'
    valuefile = out / 'plot_values.tsv.gz'
    shortlistfile = Path(__file__).with_name('boxplot_shortlist.tsv')
    hits = read_table(hitfile, ',')
    fitted = {(r['screen_id'], r['gene_id']): r for r in read_table(sourcefile)}
    values = defaultdict(list)
    for r in read_table(valuefile):
        values[(r['screen'], r['gene_id'])].append(r)
    assert set(values) == {(h['screen'], h['outcome_gene_id']) for h in hits}
    shortlist = {r['gene_symbol']: r for r in read_table(shortlistfile)}
    assert set(shortlist) <= {h['outcome'] for h in hits}
    index, selected = [], []
    with PdfPages(plots / 'candidate_boxplots_all_110.pdf') as atlas, PdfPages(plots / 'candidate_boxplots_all_110_donor_means.pdf') as donor_atlas:
        for page, hit in enumerate(hits, 1):
            key = (hit['screen'], hit['outcome_gene_id'])
            rows, fit = values[key], fitted[key]
            assert len(rows) == int(hit['n_observations'])
            assert len({r['donor'] for r in rows}) == int(hit['n_donors'])
            assert len({(r['donor'], r['SpD']) for r in rows}) == len(rows)
            for means, pdf in ((False, atlas), (True, donor_atlas)):
                fig = make_figure(hit, fit, rows, means)
                pdf.savefig(fig)
                plt.close(fig)
            entry = {'page_all': page, **hit, 'selected': int(hit['outcome'] in shortlist)}
            index.append(entry)
            if entry['selected']:
                selected.append({**entry, **shortlist[hit['outcome']]})
    selected.sort(key=lambda r: (int(r['priority']), int(r['page_all'])))
    with PdfPages(plots / 'candidate_boxplots_selected.pdf') as atlas, PdfPages(plots / 'candidate_boxplots_selected_donor_means.pdf') as donor_atlas:
        for page, hit in enumerate(selected, 1):
            key = (hit['screen'], hit['outcome_gene_id'])
            hit['page_selected'] = page
            stem = hit['screen'] + '__' + hit['outcome']
            hit['plot_file'] = 'selected/' + stem + '.pdf'
            for means, pdf in ((False, atlas), (True, donor_atlas)):
                fig = make_figure(hit, fitted[key], values[key], means, hit['rationale'])
                pdf.savefig(fig)
                if not means:
                    for extension in ('pdf', 'svg', 'png'):
                        path = plots / 'selected' / (stem + '.' + extension)
                        fig.savefig(path, dpi=300, facecolor='white')
                        if extension == 'svg':
                            path.write_text('\n'.join(line.rstrip() for line in path.read_text().splitlines()) + '\n')
                plt.close(fig)
    write_table(out / 'pair_plot_index.tsv', index, list(index[0]))
    write_table(out / 'selected_pairs.tsv', selected, list(selected[0]))
    provenance = {'inputs_sha256': {str(p.relative_to(run)): sha256(p) for p in (hitfile, sourcefile, valuefile, out / 'model_verification.tsv')},
        'script_sha256': sha256(Path(__file__)), 'export_script_sha256': sha256(Path(__file__).with_name('export_candidate_boxplots.R')),
        'shortlist_sha256': sha256(shortlistfile), 'software': {'matplotlib': matplotlib.__version__, 'numpy': np.__version__},
        'pairs': len(hits), 'selected_pairs': len(selected), 'selected_genes': len(shortlist),
        'adjustment': 'subtract centered nuisance-design contributions using saved GLS coefficients, preserving Dx; joint additionally removes M',
        'display_units': ['matched donor-SpD observations', 'within-donor means over available matched SpDs'],
        'statistics': 'unchanged saved donor-blocked limma estimates, nominal p-values and original full-family BH FDR; no boxplot tests',
        'erc_reference': 'LFF_spatial_ERC/code/22_Mediation/02_E2_E4_med_plots.Rmd: plotMedComp and cleaningY explanation',
        'outputs_sha256': {str(p.relative_to(plots)): sha256(p) for p in sorted(plots.rglob('*')) if p.is_file()}}
    (out / 'provenance.json').write_text(json.dumps(provenance, indent=2) + '\n')
    print(f'Generated {len(hits)} pair pages in each full atlas; {len(selected)} selected pairs across {len(shortlist)} genes.')


if __name__ == '__main__':
    main()
