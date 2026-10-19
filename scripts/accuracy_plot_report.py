#!/usr/bin/env python3
"""Render inspectable scientific plots from an accuracy_eval report (no refitting)."""
import argparse
import csv
import hashlib
import html
import json
from pathlib import Path

import numpy as np


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def render(report_path, output=None):
    # Import only when plotting; offline evaluation does not require matplotlib.
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    report_path = Path(report_path).resolve()
    report = json.loads(report_path.read_text(encoding='utf-8'))
    diagnostic_a = report.get('report_type') == 'mapping_a_aligned_diagnostic'
    if diagnostic_a:
        if (report.get('status') != 'A_ONLY_ALIGNED_DIAGNOSTIC'
                or report.get('independent_accuracy_assessed') is not False
                or report.get('localization_accuracy_assessed') is not False
                or any(e.get('acceptance', {}).get('pass') is not None for e in report.get('estimates', []))):
            raise ValueError('mapping A diagnostics must explicitly exclude independent accuracy and acceptance')
    elif report.get('report_type') != 'reference_disagreement':
        raise ValueError('expected a reference_disagreement or mapping_a_aligned_diagnostic report')
    out = Path(output) if output else report_path.parent / 'plots'
    if out.exists():
        raise ValueError('plot directory exists; choose a new output directory')
    inputs, datasets = {}, []
    for estimate in report['estimates']:
        for key in ('sample_csv', 'coverage_csv'):
            path = (report_path.parent / estimate[key]).resolve()
            if path.parent != report_path.parent or not path.is_file():
                raise ValueError('report CSV must be a direct file in the report directory')
            inputs[path.name] = digest(path)
        with (report_path.parent / estimate['sample_csv']).open(encoding='utf-8', newline='') as f:
            pairs = [row for row in csv.DictReader(f) if row['status'] == 'paired']
        with (report_path.parent / estimate['coverage_csv']).open(encoding='utf-8', newline='') as f:
            coverage = list(csv.DictReader(f))
        datasets.append((estimate, pairs, coverage))
    out.mkdir(parents=True, exist_ok=False)
    plt.rcParams.update({'font.size': 9, 'axes.spines.top': False, 'axes.spines.right': False,
                         'figure.dpi': 160, 'savefig.dpi': 200})
    manifest = dict(schema_version=1, report_type=report['report_type'], report_sha256=digest(report_path), input_csv_sha256=inputs,
                    interpretation='Paired samples only; no alignment refitting or extrapolation during plotting. Coverage is separate.',
                    figures=[])
    galleries = []
    all_conditions = list(dict.fromkeys(condition for estimate in report['estimates']
                                      for condition in estimate.get('conditions', [])))
    synthetic = any('synthetic' in condition.lower() for condition in all_conditions)
    sensitivity = any('LEVER_ARM_SENSITIVITY' in condition for condition in all_conditions)

    def save(fig, name):
        banners = []
        if diagnostic_a:
            banners.append('A-ONLY ALIGNED DIAGNOSTIC - NOT INDEPENDENT ACCURACY OR ACCEPTANCE')
        if synthetic:
            banners.append('SYNTHETIC TEST DATA - NOT A VEHICLE ACCURACY RESULT')
        if sensitivity:
            banners.append('LEVER-ARM SENSITIVITY - NOT AN ACCURACY ACCEPTANCE RESULT')
        elif all_conditions and not synthetic:
            banners.append('CONDITIONAL REFERENCE ASSESSMENT - SEE REPORT CONDITIONS')
        if banners:
            # Place the combined banner outside the axes/title area. A fixed
            # fraction inside short coverage figures would overlap their title.
            fig.text(.5, 1.01, '\n'.join(banners), ha='center', va='bottom',
                     color='darkred', fontsize=9,
                     bbox=dict(facecolor='white', edgecolor='darkred', alpha=.95))
        fig.savefig(out / (name + '.png'), bbox_inches='tight')
        fig.savefig(out / (name + '.pdf'), bbox_inches='tight')
        plt.close(fig)
        manifest['figures'].append(dict(name=name, png_sha256=digest(out / (name+'.png')),
                                       pdf_sha256=digest(out / (name+'.pdf'))))
        return name + '.png'

    for number, (estimate, pairs, coverage) in enumerate(datasets):
        name = estimate['name']
        prefix = f'{number+1:02d}'
        images = []
        fig, ax = plt.subplots(figsize=(7, 5))
        if pairs:
            column = lambda key: np.array([float(row[key]) for row in pairs])
            time_ns = [int(row['timestamp_ns']) for row in pairs]
            t = np.array([(stamp-time_ns[0])/1e9 for stamp in time_ns])
            ax.scatter(column('reference_x'), column('reference_y'), s=5, alpha=.6, label='CGI reference')
            ax.scatter(column('estimate_enu_x'), column('estimate_enu_y'), s=4, alpha=.6, label=name)
            ax.legend()
        else:
            ax.text(.5, .5, 'No qualified paired samples', ha='center', transform=ax.transAxes)
        ax.set(xlabel='East [m]', ylabel='North [m]', title=f'{name}: fixed ENU, rear axle')
        ax.set_aspect('equal', adjustable='datalim')
        ax.grid(alpha=.25)
        images.append(save(fig, prefix+'_trajectory_xy'))
        for suffix, keys, units in (
                ('position', ['x_m', 'y_m', 'z_m', 'horizontal_m', 'position_3d_m', 'lateral_m'], 'm'),
                ('attitude', ['roll_deg', 'pitch_deg', 'yaw_deg', 'attitude_geodesic_deg'], 'deg')):
            fig, axes = plt.subplots(len(keys), 1, figsize=(9, 1.6*len(keys)), sharex=True)
            labels = {'x_m': 'X error', 'y_m': 'Y error', 'z_m': 'Z error',
                      'horizontal_m': 'Horizontal', 'position_3d_m': '3D position', 'lateral_m': 'Lateral FLU',
                      'roll_deg': 'Roll', 'pitch_deg': 'Pitch', 'yaw_deg': 'Yaw',
                      'attitude_geodesic_deg': 'SO(3) angle'}
            for ax, key in zip(axes, keys):
                if pairs:
                    values = column(key)
                    ax.scatter(t, values, s=4, alpha=.65)
                    if np.ptp(values) < 1e-10:
                        margin = .01 if units == 'm' else .05
                        ax.set_ylim(float(np.mean(values))-margin, float(np.mean(values))+margin)
                ax.ticklabel_format(style='plain', axis='y', useOffset=False)
                ax.set_ylabel(f'{labels[key]} [{units}]')
                ax.grid(alpha=.25)
            axes[0].set_title(f'{name}: estimate minus reference; paired samples only')
            axes[-1].set_xlabel('Seconds from first paired sample (inspect coverage for outages)')
            images.append(save(fig, prefix+'_'+suffix))
        fig, ax = plt.subplots(figsize=(10, 3.2))
        states = sorted({row['state'] for row in coverage})
        epochs = sorted({row['epoch'] for row in coverage})
        palette = plt.get_cmap('tab10')
        colors = {state: palette(i % 10) for i, state in enumerate(states)}
        for epoch_index, epoch in enumerate(epochs):
            segments = [row for row in coverage if row['epoch'] == epoch]
            start = min(int(row['start_ns']) for row in segments)
            labeled = set()
            for row in segments:
                state = row['state']
                ax.broken_barh([((int(row['start_ns'])-start)/1e9, float(row['duration_sec']))],
                               (epoch_index-.35, .7), facecolors=colors[state],
                               label=state if state not in labeled and epoch_index == 0 else None)
                labeled.add(state)
        ax.set(xlabel='Seconds from declared interval start in each epoch', ylabel='Clock epoch',
               title=f'{name}: full recording coverage (not a pose-error plot)',
               yticks=range(len(epochs)), yticklabels=epochs)
        if states:
            ax.legend(loc='upper center', bbox_to_anchor=(.5, -.2), ncol=min(4, len(states)))
        ax.grid(axis='x', alpha=.25)
        images.append(save(fig, prefix+'_coverage'))
        galleries.append((name, estimate['trajectory_kind'], estimate['acceptance']['status'], images))
    with (out / 'plot_manifest.json').open('x', encoding='utf-8') as f:
        json.dump(manifest, f, indent=2, allow_nan=False)
    content = ['<!doctype html><meta charset="utf-8"><title>Trajectory reference assessment</title>',
               '<style>body{max-width:1100px;margin:30px auto;font:16px system-ui;line-height:1.6}img{max-width:100%}section{margin:3em 0}</style>',
               '<h1>Trajectory reference assessment</h1><p>Fixed alignment; CGI-430 is a quality-gated reference. '
               'These figures show reference disagreement, with recording coverage reported separately. '
               'SO(3) geodesic is the main attitude metric; Euler differences are diagnostics.</p>']
    if synthetic:
        content.append('<p style="color:darkred;font-weight:bold">SYNTHETIC TEST DATA — NOT A VEHICLE ACCURACY RESULT.</p>')
    if diagnostic_a:
        content.append('<p style="color:darkred;font-weight:bold">A-ONLY ALIGNED DIAGNOSTIC — '
                       'alignment was fitted using the same mapping data. '
                       'These are conditional residuals, not independent absolute accuracy, '
                       'map localization accuracy, or acceptance results.</p>')
    if sensitivity:
        content.append('<p style="color:darkred;font-weight:bold">LEVER-ARM SENSITIVITY — NOT AN ACCURACY ACCEPTANCE RESULT. No best lever arm is selected.</p>')
    if all_conditions:
        content.append('<h2>Conditions and evidence limits</h2><ul>')
        content.extend('<li>' + html.escape(condition) + '</li>' for condition in all_conditions)
        content.append('</ul>')
    for name, kind, status, images in galleries:
        content.append(f'<section><h2>{html.escape(name)}</h2><p>{html.escape(kind)} / {html.escape(status)}</p>')
        content += [f'<p><a href="{img[:-4]}.pdf">PDF</a></p><img alt="{html.escape(img)}" src="{img}">' for img in images]
        content.append('</section>')
    content.append('<p><a href="plot_manifest.json">Plot and source hashes</a></p>')
    (out / 'index.html').write_text('\n'.join(content), encoding='utf-8')
    return dict(output=str(out), figures=len(manifest['figures']))


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--report', required=True, help='accuracy_eval report.json')
    p.add_argument('--output', help='new directory; defaults to report/plots')
    a = p.parse_args()
    try:
        print(json.dumps(render(a.report, a.output)))
    except (ValueError, KeyError, OSError) as error:
        p.exit(2, f'accuracy_plot_report: {error}\n')


if __name__ == '__main__':
    main()
