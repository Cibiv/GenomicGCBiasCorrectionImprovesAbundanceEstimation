#!/usr/bin/env python3
# %%
import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import os
import sys
from matplotlib.ticker import FuncFormatter
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.pipeline import make_pipeline

plt.rcParams['font.sans-serif'] = ['Helvetica']
plt.rcParams['font.family'] = 'sans-serif' # Set the generic family to 'sans-serif' to use the specified font


# Set the global font size
plt.rcParams['font.size'] = 8

def cm_to_inches(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)

path = 'guacamole_results/'

print("Initializing", file=sys.stderr)
sra_all = pd.read_csv("sra_ids.txt", sep='\t', header=None, usecols=[0, 3])
sra_all.columns = ['sra_id', 'sample']

files = os.listdir(path)
mpa_files = os.listdir('metaphlan_results')

samples = [file.split('_')[2] for file in files]

true_abundances = [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 2, 1, 1, 1, 1]
true_abundances = true_abundances / np.sum(true_abundances)
taxids = [303, 562, 820, 823, 907, 1282, 1309,
          1423, 1531, 1584, 1747, 33035, 39491, 74426, 216816, 239935,
          292800, 411483, 437897]
true_ab = dict(zip(taxids, true_abundances))

tmp_df = pd.read_csv(path+'output_A0_1.txt', sep='\t', index_col=False)
df = tmp_df.set_index('taxonomy_id', drop=True)

true_species = [df.loc[id]['name'] for id in taxids]

# Load Bracken/GuaCAMOLE output to map TaxIDs to species names
bracken_file = os.path.join(path, "output_A0_1.txt")
bracken_df = pd.read_csv(bracken_file, sep="\t")

# Create a dictionary mapping species names to NCBI TaxIDs
species_to_taxid = dict(zip(bracken_df["name"], bracken_df["taxonomy_id"]))

# Fix GTDB → NCBI mismatch
species_to_taxid["Blautia coccoides"] = 33035
species_to_taxid["Pseudomonas_E putida"] = 303
species_to_taxid["Escherichia coli_F"] = 562
species_to_taxid["Escherichia sp005843885"] = 562
species_to_taxid["Bifidobacterium infantis"] = 216816
species_to_taxid["Lactobacillus sp027681045"] = 1584
species_to_taxid["Lactobacillus sp944327175"] = 1584
# motus
species_to_taxid["[Eubacterium] rectale"] = 39491
species_to_taxid["[Clostridium] clostridioforme/bolteae"] = 1531
species_to_taxid["Megamonas funiformis/rupellensis"] = 437897

# Initialize lists for dataframe creation
error = []
samples = []
methods = []
txid_error = []
method_txid = []
sample_txid = []
abundance = []
txid = []
f_positive = []
txid_fp = []
fp_samples = []
fp = []
abundance_fp = []
replicates = []
replicate_txid = []
frac_fp = []
removed = []

# Define methods and their corresponding column names
methods_dict = {
    "GuaCAMOLE": "GuaCAMOLE_estimate",
    "Bracken": "Bracken_estimate",
    "MetaPhlAn4": "mpa",  # Special handling for MetaPhlAn4
    "Sylph": "Taxonomic_abundance",  # Extract abundances from Sylph output
    "SingleM": "taxonomy",  # Needs parsing
    "mOTUS": "species"  # Needs species-to-TaxID conversion
}
replicate = 1
prev_sample = None

print("Reading results of for individual tools", file=sys.stderr)
for nrow in range(sra_all.shape[0]):
    sample = sra_all['sample'][nrow]
    sra = sra_all['sra_id'][nrow]
    print(f"  {sample}: {sra}", file=sys.stderr)

    if prev_sample == sample:
        replicate += 1
    else:
        replicate = 1

    prev_sample = sample

    if not os.path.exists(os.path.join(path, 'output_' + sample + '_' + str(replicate) + '.txt')):
        continue

    # Load MetaPhlAn4 results
    mpa_path = f'metaphlan_results/{sample[0]}_{sample[1]}_{replicate}.mpa'
    mpa_df = pd.read_csv(mpa_path, sep='\t', comment='#', header=None)
    mpa_df = mpa_df.loc[mpa_df[0].str.contains('s__') & ~mpa_df[0].str.contains('t__')]
    mpa_txids = np.array([int(i[6]) for i in mpa_df[1].str.split('|')])
    mpa_abundances = dict(zip(mpa_txids, np.array(mpa_df[2])))
    mpa_abundances[411483] = mpa_abundances[853]
    del mpa_abundances[853]

    # Load GuaCAMOLE and Bracken estimates
    tmp_df = pd.read_csv(f'{path}output_{sample}_{replicate}.txt', sep='\t', index_col=False)
    df = tmp_df.set_index('taxonomy_id', drop=True)

    ab_true = np.array([true_ab[taxid] for taxid in np.sort(list(true_ab.keys()))])

    f_positive.extend(['true_positive' if taxid in true_ab.keys() else 'false_positive' for taxid in df.index])
    txid_fp.extend(df.index)
    abundance_fp.extend(df['Bracken_estimate'])
    frac_fp.append(np.sum(df.loc[df.index.isin(true_ab.keys()), 'Bracken_estimate']))
    removed.append(np.sum(df.loc[df['GuaCAMOLE_estimate'] == 0, 'Bracken_estimate']))
    fp_samples.extend(np.repeat(sample, df.shape[0]))

    for method_name, column_name in methods_dict.items():
        if method_name == "MetaPhlAn4":
            ab_method = np.array(
                [mpa_abundances[taxid] if taxid in mpa_abundances.keys() else 0 for taxid in true_ab.keys()]
            ) / 100

        elif method_name in ["GuaCAMOLE", "Bracken"]:
            ab_method = np.array(
                [df[column_name][taxid] if taxid in df.index else 0 for taxid in true_ab.keys()]
            )
            ab_method /= np.sum(ab_method)

        elif method_name == "Sylph":
            sylph_df = pd.read_csv(f"sylph_results/{sra}_{sample[0]}_{sample[1]}_{replicate}", sep="\t")
            sylph_df['species'] = [' '.join(name[1:3]) for name in sylph_df['Contig_name'].str.split(' ')]
            sylph_df["NCBI_taxid"] = sylph_df["species"].apply(lambda x: species_to_taxid.get(x.split("/")[-1], None))
            sylph_df = sylph_df.groupby("NCBI_taxid", as_index=False).agg({
                "Taxonomic_abundance": "sum",
                "NCBI_taxid": "first",
                "species": "first"  # Keeps the first occurrence
            })
            sylph_df = sylph_df.dropna(subset=["NCBI_taxid"])  # Remove unmatched taxa

            ab_method = np.array(
                [sylph_df.loc[sylph_df["NCBI_taxid"] == taxid, "Taxonomic_abundance"].sum()
                 if taxid in true_ab.keys() else 0 for taxid in true_ab.keys()]
            )
            ab_method /= np.sum(ab_method)

        elif method_name == "SingleM":
            singlem_df = pd.read_csv(f"singlem_results/{sra}_{sample[0]}_{sample[1]}_{replicate}_profile.txt",
                                     sep="\t", skiprows=[1])
            singlem_df['species'] = [name.split("; ")[-1].replace("s__", "") for name in singlem_df['taxonomy']]
            singlem_df.columns = ['taxonomy', 'abundance', 'species']
            singlem_df['abundance'] = singlem_df['abundance'].astype(float)

            singlem_df["NCBI_taxid"] = singlem_df["species"].apply(lambda x: species_to_taxid.get(x, None))
            singlem_df = singlem_df.dropna(subset=["NCBI_taxid"])
            singlem_df['NCBI_taxid'] = singlem_df['NCBI_taxid'].astype(int)
            singlem_df = singlem_df.groupby("NCBI_taxid", as_index=False).agg({
                "abundance": "sum",
                "NCBI_taxid": "first",
                "species": "first"  # Keeps the first occurrence
            })

            ab_method = np.array(
                [singlem_df.loc[singlem_df["NCBI_taxid"] == taxid, "abundance"].sum()
                 if taxid in true_ab.keys() else 0 for taxid in true_ab.keys()]
            )
            ab_method /= np.sum(ab_method)

        elif method_name == "mOTUS":
            motus_df = pd.read_csv(f"motus_results/{sra}_{sample[0]}_{sample[1]}_{replicate}.bz2", sep="\t", skiprows=2)
            motus_df.columns = ["species_full", "abundance"]

            # Extract only the first two words of species name
            motus_df["species"] = motus_df["species_full"].str.split().str[:2].str.join(" ")

            # Remove zero-abundance species
            motus_df = motus_df[motus_df["abundance"] > 0]

            # Map species to NCBI TaxIDs
            motus_df["NCBI_taxid"] = motus_df["species"].apply(lambda x: species_to_taxid.get(x, None))
            motus_df = motus_df.dropna(subset=["NCBI_taxid"])
            motus_df["NCBI_taxid"] = motus_df["NCBI_taxid"].astype(int)


            ab_method = np.array(
                [motus_df.loc[motus_df["NCBI_taxid"] == taxid, "abundance"].sum()
                 if taxid in true_ab.keys() else 0 for taxid in true_ab.keys()]
            )
            ab_method /= np.sum(ab_method)

        re_method = (ab_method - ab_true) / ab_true
        #lf_method = np.log(ab_method / ab_true)

        txid_error.extend(re_method)
        abundance.extend(ab_method)
        method_txid.extend(np.repeat(method_name, len(re_method)))
        sample_txid.extend(np.repeat(sample, len(re_method)))
        replicate_txid.extend(np.repeat(replicate, len(re_method)))
        error.append(np.mean(abs(re_method)))
        methods.append(method_name)
        samples.append(sample)
        replicates.append(replicate)
        txid.extend(list(true_ab.keys()))


# %%
err_df = pd.DataFrame({
    'error': error, 'method': methods, 'sample': samples, 'replicate': replicates
})

txid_df = pd.DataFrame({
    'abundance': abundance, 'method': method_txid, 'sample': sample_txid, 'txid': txid, 'error': txid_error, 'replicate': replicate_txid
})

order_libprep = ['A0', 'AL', 'AH', 'B0', 'BL', 'BH', 'C0', 'CL', 'CH',
                 'D0', 'DL', 'DH', 'E0', 'EL', 'EH', 'F0', 'FL', 'FH',
                 'GL', 'GH', 'HL', 'HH', 'IL', 'IH', 'JL', 'JH', 'KL',
                 'KH']

markers_dict = {
    'GuaCAMOLE': 'X',   # Large X (stays the same)
    'Bracken': 'o',     # Circle (stays the same)
    'MetaPhlAn4': 's',  # Square (stays the same)
    'SingleM': 'P',     # Valid alternative → Filled Plus (bold +)
    'mOTUS': 'D',       # Valid alternative → Diamond
    'Sylph': 'v'        # Valid alternative → Triangle Down
}
hue_order_list = ['GuaCAMOLE', 'Bracken', 'MetaPhlAn4', 'SingleM', 'Sylph', 'mOTUS']

err_df['Method'] = err_df['method']
err_df['error'] = err_df['error'] * 100
err_df['sample'] = pd.Categorical(err_df['sample'], categories=order_libprep, ordered=True)

# Sort the DataFrame by the new categorical 'sample' column
err_df = err_df.sort_values('sample')
# err_df = err_df.loc[(err_df['method'] != 'GuaCAMOLE'), :]
plt.figure(figsize=cm_to_inches(12, 6))
g = err_df.groupby(['sample', 'Method'], observed=False)
err_df_avg = g.mean('error').reset_index()
jitter_strength = 0.1
err_df['sample_jitter'] = err_df['sample'].cat.codes + np.random.uniform(-jitter_strength, jitter_strength, size=len(err_df))
ax = sns.scatterplot(x="sample_jitter", y="error", data=err_df, style='Method', markers=markers_dict,
                     s=10, hue='Method', hue_order=hue_order_list)

plt.vlines([2.5, 5.5, 8.5, 11.5, 14.5, 17.5, 19.5, 21.5, 23.5, 25.5], ymin=0,
           ymax=np.max(err_df['error']), alpha=0.5,
           linewidth=0.5, colors='grey')

plt.xticks(ticks=range(len(order_libprep)), labels=order_libprep)
plt.gca().set_xticklabels([label.get_text()[1] for label in plt.gca().get_xticklabels()])
plt.xlabel('Library Preparation Protocol')
plt.ylabel('Relative Error (%)')
plt.legend(loc='upper left')
plt.ylim(0, max(err_df['error']) + 1)
plt.tight_layout()
plt.savefig('avg_error_vs_protocol.pdf', transparent=True)
plt.close()

# %% boxplot with lines
plt.figure(figsize=cm_to_inches(6, 6))

jitter = 0.05
df_x_jitter = np.random.normal(loc=0, scale=jitter, size=err_df_avg.shape[0])
method_int = np.zeros(len(err_df_avg['error']))
method_int[err_df_avg['Method'] == 'Bracken'] = 0
method_int[err_df_avg['Method'] == 'GuaCAMOLE'] = 1
method_int[err_df_avg['Method'] == 'MetaPhlAn4'] = 2
method_int[err_df_avg['Method'] == 'SingleM'] = 3
method_int[err_df_avg['Method'] == 'Sylph'] = 4
method_int[err_df_avg['Method'] == 'mOTUS'] = 5
method_int += df_x_jitter
err_df_avg['method_int'] = method_int

fig, ax = plt.subplots(figsize=cm_to_inches(6, 6))
PROPS = {
    'boxprops': {'facecolor': 'none', 'edgecolor': 'black'},
    'medianprops': {'color': 'black'},
    'whiskerprops': {'color': 'black'},
    'capprops': {'color': 'black'}
}

marked_width = 2

sample_dict = {'GH': (sns.color_palette()[3], marked_width),
               'DH': (sns.color_palette()[4], marked_width),
               'IH': (sns.color_palette()[5], marked_width),
               'FH': (sns.color_palette()[9], marked_width),
               'IL': (sns.color_palette()[8], marked_width)}

sns.boxplot(x='Method', y='error', data=err_df_avg, showfliers=False, **PROPS)
for method in np.unique(err_df_avg['Method']):
    for sample in np.unique(err_df_avg.loc[err_df_avg['Method'] == method, 'sample']):
        if sample in sample_dict:
            color = sample_dict[sample][0]
        else:
            color = 'grey'
        ind = (err_df_avg['Method'] == method) & (err_df_avg['sample'] == sample)
        ax.plot(err_df_avg.loc[ind, 'method_int'],
                err_df_avg.loc[ind, 'error'], 'o', ms=2, color=color)


plt.ylim(0, max(err_df_avg['error']) + 1)
plt.xticks(rotation=45)
plt.ylabel('Relative Error (%)')
plt.xlabel('')
plt.tight_layout()
plt.savefig('boxplot.pdf', transparent=True)
plt.close()

# %% GC content vs Error with quadratic regression
ref_dist = np.loadtxt('ref_bin_100_input.dist')
txid_df_avg = txid_df.groupby(['method', 'sample', 'txid']).mean('error').reset_index()
gc_content = np.argmax(ref_dist, axis=0) + np.random.uniform(0.1, 0.02, size=ref_dist.shape[1])
gc_dct = dict(zip(np.sort(taxids).tolist(), gc_content.tolist()))
txid_df_avg['GC content'] = [float(gc_dct[taxid]) for taxid in txid_df_avg['txid']]
txid_df_avg['PCR'] = [sample[1] for sample in txid_df_avg['sample']]
txid_df['GC content'] = [float(gc_dct[taxid]) for taxid in txid_df['txid']]


methods = txid_df['method'].unique()
subset = txid_df.loc[txid_df['sample'] == 'GH'].copy()
subset['error'] = abs(subset['error']) * 100

g = sns.lmplot(x='GC content', y='error', hue='method', data=subset, legend=False,
           order=2, scatter=False, ci=None, height=cm_to_inches(6, 6)[0], aspect=11/10)

for ax in g.axes.flat:
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

#plt.ylim(0, 30)
plt.ylabel('Relative Error (%)')
plt.xlabel('GC content (%)')
sns.scatterplot(x='GC content', y='error', hue='method', data=subset, s=10, style='method',
                legend=True, hue_order=hue_order_list, markers=markers_dict)
plt.tight_layout()
plt.savefig('regplot_GH_quad.pdf', transparent=True)
plt.close()
# %% PCR Cycles vs Error with regression

pcr_cycles_dict = {
    'A0': 0, 'AL': 4, 'AH': 9,
    'B0': 0, 'BL': 8, 'BH': 12,
    'C0': 0, 'CL': 8, 'CH': 8,
    'D0': 0, 'DL': 5, 'DH': 15,
    'E0': 0, 'EL': 4, 'EH': 14,
    'F0': 0, 'FL': 4, 'FH': 14,
    'G0': None, 'GL': 5, 'GH': 15,
    'H0': None, 'HL': 4, 'HH': 9,
    'I0': None, 'IL': 4, 'IH': 9,
    'J0': None, 'JL': 4, 'JH': 12,
    'K0': None, 'KL': 6, 'KH': 11
}
err_df['cycles'] = [pcr_cycles_dict[sample] for sample in err_df['sample']]
err_df['abs_error'] = abs(err_df['error'])

g = sns.lmplot(x='cycles', y='abs_error', hue='method', data=err_df, order=2, hue_order = hue_order_list,
           scatter=False, ci=None, legend=False, height=cm_to_inches(6, 6)[0], aspect=11/10)

for ax in g.axes.flat:
    ax.spines['top'].set_visible(True)
    ax.spines['right'].set_visible(True)
    ax.spines['left'].set_visible(True)
    ax.spines['bottom'].set_visible(True)

plt.gca().xaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{int(x)}'))
plt.ylabel('Relative Error (%)')
plt.xlabel('PCR cycles')
sns.scatterplot(x='cycles', y='abs_error', hue='method', data=err_df, s=10, style='method',
                legend=True, hue_order=hue_order_list, markers=markers_dict)
plt.tight_layout()
plt.savefig('regplot_pcr.pdf', transparent=True)
plt.close()
