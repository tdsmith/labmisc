#!/usr/bin/env python
# coding=utf8

# spec_agg.py
# Aggregates wave scans produced by the George Lab's spectrophotometer.
# Reads all .csv files in the current directory and creates a series of .png
# charts and spec_summary.html. Does not check before overwriting files so be
# careful. Also spits a table of A230, A260, A280, and A320 as tab-separated
# values to stdout. A convenient way to print the graphs for your lab notebook
# is to open the .html file in Word, select Print..., and set scaling to 40%.
# Invoke as: ./spec_agg.py > output.tab

from glob import glob
import codecs
import matplotlib as mpl
from matplotlib import pyplot as plt
from os.path import splitext

def fuzzy_lookup(needle, haystack, tol=0.01):
    # naïve O(N) method to look up a float value
    for i in haystack:
        if abs(needle - i[0]) < tol: return i[1]
    raise ValueError

df = glob('*.csv')

d = {}
max_abs = 0

for fn in df:
    f = codecs.open(fn, 'r', 'utf-16') # your guess is as good as mine
    buf = f.readlines()[10:] # trim header
    buf = [line[:-1].split(',') for line in buf] # ignore trailing \n
    scan = [(float(line[0]), float(line[1])) for line in buf if len(line) > 1]
    max_abs = max(max_abs, max([row[1] for row in scan]))
    d[splitext(fn)[0]] = scan

html = open('spec_summary.html', 'w')
print >> html, """<!DOCTYPE html>
<html>
<head><title>Spectrophotometer output</title></head>
<body>
<table><tbody>
<tr>""",

for i, k in enumerate(sorted(d.keys())):
    if i % 3 == 0 and i > 0: print >> html, '</tr>\n<tr>'
    fig = plt.figure(figsize=(4,2),
            subplotpars=mpl.figure.SubplotParams(left=0.18, bottom=0.25, right=0.95, top=0.85))
    ax = fig.gca()
    ax.plot(*zip(*d[k]))
    ax.set_ylim(0, max_abs)
    ax.set_title(k)
    ax.set_xlabel('Wavelength (nm)')
    ax.set_ylabel('Absorbance')
    fig.savefig('%s.png' % k)
    print >> html, '<td><img src="%s.png" /></td>' % k,

print >> html, """</tr>
</tbody></table>
</body>
</html>"""

html.close()

print 'Label\tA230\tA260\tA280\tA320'
for k in sorted(d.keys()):
    e = dict(zip(['a','b','c','d'], [fuzzy_lookup(needle, d[k]) for needle in [230, 260, 280, 320]]))
    e['k'] = k
    print '%(k)s\t%(a).3f\t%(b).3f\t%(c).3f\t%(d).3f' % e


