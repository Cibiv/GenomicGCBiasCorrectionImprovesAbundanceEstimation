#!/usr/bin/env python3
import numpy as np
import pandas as pd
import os
import fnmatch as fnm
import seaborn as sns
from matplotlib import pyplot as plt

plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['font.family'] = 'sans-serif' # Set the generic family to 'sans-serif' to use the specified font

def cm_to_inches(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)

eff_all = [f for f in os.listdir('guacamole_results') if fnm.fnmatch(f, "efficiencies_*.txt") ]

samples = []
replicates = []
gc_content = []
efficiencies = []
for file in eff_all:
    eff = np.loadtxt(os.path.join('guacamole_results', file))
    sample = file[-8:-6]
    replicate = file[-5]
    samples.extend(np.repeat(sample, len(eff)))
    replicates.extend(np.repeat(replicate, len(eff)))
    gc_content.extend(np.arange(len(eff)))
    eff = eff / max(eff)
    efficiencies.extend(eff)

df = pd.DataFrame(
    {
        'Efficiency': efficiencies,
        'Sample': samples,
        'GC content': gc_content,
        'Replicate': replicates
    }
)
df = df.loc[df['Efficiency'] != 0, :]

marked_width = 2

sample_dict = {'GH': (sns.color_palette()[3], marked_width),
               'DH': (sns.color_palette()[4], marked_width),
               'IH': (sns.color_palette()[5], marked_width),
               'FH': (sns.color_palette()[9], marked_width),
               'IL': (sns.color_palette()[8], marked_width)}

palette = sns.color_palette(['grey'], len(np.unique(df['Sample'])))
plt.figure(figsize=cm_to_inches((6, 6)))
plt.rcParams.update({'font.size': 12})
ax = sns.lineplot(x='GC content', y='Efficiency', hue='Sample', palette=palette, data=df, errorbar=None, linewidth=0.5)
plt.legend([], [], frameon=False)
plt.ylim(0, 1.1)
plt.ylabel("Sequencing Efficiency")
plt.xlabel("GC content (%)")



for sample in sample_dict:
    g = df.loc[df['Sample'] == sample, 'Efficiency'].groupby(df.loc[df['Sample'] == sample, 'GC content'])
    eff = g.mean()
    ax.plot(eff.index, eff, color=sample_dict[sample][0], linewidth=sample_dict[sample][1])

plt.savefig('efficiencies.pdf', bbox_inches="tight")
