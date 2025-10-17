# %%
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import os
import pandas as pd


def cm_to_inches(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)


# should be one of ['increasing', 'decreasing', 'concave']
bias = 'decreasing'


def efficiency(gc, bias):
    if bias == "increasing":
        return 1.25 * gc
    elif bias == "decreasing":
        return -1.25 * gc + 1
    elif bias == "concave":
        return 1 - 10 * (gc - 0.5) ** 2


eff_theory = [efficiency(gc=x / 100, bias=bias) for x in range(100)]
eff_method = np.loadtxt(bias + '/efficiencies.txt')
meth_max = np.argmax(eff_method)
eff_theory = np.array(eff_theory) / eff_theory[meth_max]
eff_method = eff_method / max(eff_method)

plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(figsize=cm_to_inches(6, 6))
sns.lineplot(x=range(100), y=eff_theory, color='red', linestyle='--', linewidth=3, label='True Efficiency')
sns.lineplot(x=range(101), y=eff_method, alpha=.7, linewidth=3, color='blue', label='GuaCAMOLE')
plt.xlim(26, 68)
plt.ylim(0, 1.1)
plt.xlabel('GC content (%)')
plt.ylabel('Relative Sequencing Efficiency')
plt.legend(frameon=False)
plt.savefig('fig2a/' + bias + '_eff.pdf')

## barplots

# %%
true_abundances = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]
true_abundances = true_abundances / np.sum(true_abundances)
taxids = [303, 562, 820, 823, 907, 1282, 1309,
          1423, 1531, 1584, 1747, 33035, 39491, 74426, 216816, 239935,
          292800, 411483, 437897]
true_ab = dict(zip(taxids, true_abundances))

tmp_df = pd.read_csv(bias + '/res_df.csv', usecols=[1, 2, 3, 4])
df = tmp_df.set_index('taxid.1', drop=True)
df['true'] = true_ab

ref_dist = np.loadtxt('ref_bin_100_input.dist')
gc = np.argmax(ref_dist, axis=0) + np.random.uniform(0, 0.001, size=ref_dist.shape[1])

er_guac = (df['abundance_ls'] - df['true'])*100 / df['true']
er_br = (df['abundance_br'] - df['true'])*100 / df['true']

sorted_indices = np.argsort(gc)
gc_sorted = np.array(gc)[sorted_indices]
er_br_sorted = np.array(er_br)[sorted_indices]
er_guac_sorted = np.array(er_guac)[sorted_indices]

plt.rcParams.update({'font.size': 11})
fig, ax = plt.subplots(figsize=cm_to_inches(6, 6))

# Define bar width and positions
bar_positions = gc_sorted

df = pd.DataFrame(
    {
    'error': np.concatenate([er_guac_sorted, er_br_sorted]),
    'GC content': np.tile(bar_positions, 2),
    'method': np.repeat(['GuaCAMOLE', 'Bracken'], len(er_br_sorted))
    }
)
style_dict = {'GuaCAMOLE': 'X', 'Bracken': 'o'}

sns.lineplot(x='GC content', y='error', hue='method', data=df, legend=False)
sns.scatterplot(x='GC content', y='error', hue='method', data=df, style='method',markers=style_dict)

#plt.legend(frameon=False)
plt.ylim(-40, 65)
plt.xticks(rotation=45)
plt.xlabel('GC content (%)')
plt.ylabel('Relative Error (%)')
plt.xlim(30, 65)
plt.savefig('fig2a/' + bias + '_scatter.pdf')

