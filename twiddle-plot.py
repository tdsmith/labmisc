#!/usr/bin/env python
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
import sys

df = pd.read_table(sys.argv[1], delim_whitespace=True, names=["time", "n", "Vcc", "Vr", "R", "T"])
df['time'] = pd.to_datetime(df['time'], unit='s', utc=True).map(lambda x: x.tz_localize('UTC').tz_convert('America/Los_Angeles'))
p = df.plot('time', 'T', ylim=(20,40))
p.axhline(37., color='gray')
plt.show()

