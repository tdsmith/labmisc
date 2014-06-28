#!/usr/bin/env python3

from __future__ import division
import argparse
import scipy as sp
import scipy.optimize
import scipy.stats
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import sys
import os
import time

def standardize(infile, max_q, factor):
    "Fit OD = (A-D)/(1+(x/C)^B) + D"
    plate = np.loadtxt(infile, delimiter=',')
    od = np.ravel(plate[0:8, 0:2], order='F')
    conc = ([max_q*factor**(-i) for i in range(7)] + [0]) * 2
    def f(x, a, b, c, d):
            return (a-d)/(1+(x/c)**b) + d
    guess = [plate[7,0], 1, max_q/factor, plate[1,0]]
    popt, pcov = sp.optimize.curve_fit(f, conc, od, guess)

    x = np.linspace(0, max_q, 200)
    fig, ax = plt.subplots()
    ax.plot(x, f(x, *popt), '-', conc, od, 'o')
    ax.set_xlim(0, max_q*1.1)
    ax.set_xlabel('Concentration (ng/ml)')
    ax.set_ylabel('OD')
    ax.set_title(infile.name)

    def inv(od, a, b, c, d):
        return c * ((a-d)/(od-d) - 1)**(1/b)
    std_plate = inv(plate, *popt)

    pearson_r, p_val = sp.stats.pearsonr(f(conc, *popt), od)

    return std_plate, {'params': popt, 'r': pearson_r, 'p': p_val}, fig

def main():
    description=('Standardize an ELISA plate.\n'
        'The standard curve is assumed to be in columns 1 and 2. '
        'H1 and H2 are assumed to be blanks.\n'
        'Prints a standardized plate to stdout. Optionally saves summary '
        'statistics to an output file.')
    parser = argparse.ArgumentParser(description=description,
            formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--max', type=float, default=500.0,
            help='The top of the standard curve')
    parser.add_argument('--factor', type=float, default=2,
            help='Serial dilution factor for the standard curve')
    parser.add_argument('--report', help='Directory to plate reports in. '
            'Default is to skip reports. Overwrites if present.', required=False)
    parser.add_argument('input_file', help='CSV file to operate on', type=argparse.FileType('r'))
    args = parser.parse_args()

    plate, stats, graph = standardize(args.input_file, args.max, args.factor)
    np.savetxt(sys.stdout.buffer, plate, fmt='%.3f', delimiter=',')

    if not args.report:
        return

    # generate the report
    try:
        os.mkdir(args.report)
    except OSError:
        pass

    graph.savefig(os.path.join(args.report, 'standard_curve.png'))
    html = ('<html><body><p>Report for {filename}, generated at {time}'
            '<table><tr><td><img src="standard_curve.png" /></td>'
            '<td><p>OD = ({params[0]:.3g}-{params[3]:.3g})/(1 + (Concentration/'
            '{params[2]:.3g})^{params[1]:.3g}) + {params[3]:.3g}</p>'
            '<p>R<sup>2</sup>: {r2:g}</p><p><i>p</i>: {p:g}</p>'
            '</td></tr></table></body></html>').format(
                    filename=os.path.basename(args.input_file.name),
                    time=time.ctime(),
                    params=stats['params'],
                    r2=stats['r']**2,
                    p=stats['p'])
    with open(os.path.join(args.report, 'index.html'), 'w') as f:
        f.write(html)


if __name__ == '__main__':
    main()

