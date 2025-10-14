#!/usr/bin/env python3

import pandas as pd
import numpy as np
import seaborn as sns
from matplotlib import pyplot as plt
import os
import warnings # Import warnings module

cb_palette = sns.color_palette('colorblind')

method_colors = {
    "GuaCAMOLE": cb_palette[0],  # Blue
    "Bracken": cb_palette[1],  # Orange
    "MetaPhlAn4": cb_palette[2],  # Green
    "GuaCAMOLE_eff": cb_palette[3]  # Red (Example for the efficiency-based estimate)
    # Add other methods if needed
}

# Choose markers (standard matplotlib markers)
method_markers = {
    "GuaCAMOLE": "o",  # Circle
    "Bracken": "s",  # Square
    "MetaPhlAn4": "X",  # X
    "GuaCAMOLE_eff": "^"  # Triangle up
    # Add other methods if needed
}

parameters = [param.split('.')[0] for param in os.listdir('parameters')]

def cm_to_inches(*tupl):
    inch = 2.54
    if isinstance(tupl[0], tuple):
        return tuple(i / inch for i in tupl[0])
    else:
        return tuple(i / inch for i in tupl)

assemblies = pd.read_csv('assembly_summary.txt.bz2',sep='\t', header=0, skiprows=1)
assemblies.set_index('#assembly_accession', inplace=True)

# files = os.listdir('results') # No longer needed if iterating through parameters
# samples = [file.split('_')[2] for file in files] # No longer needed

# Initialize lists for dataframe creation
error = []
samples_list = []
methods = []
replicates = []
f_positives = []
f_negatives = []

txid_error = []
method_txid = []
sample_txid = []
replicate_txid = []
abundance_list = []
txid_list = []

gc_medium = []
gc_extreme = []
gc_uniform = []

efficiencies = [] # This list wasn't used later, consider removing if still unused
eff_sample = [] # This list wasn't used later, consider removing if still unused

# --- Keep track of skipped replicates ---
skipped_replicates = []
# --------------------------------------

# Define methods and their corresponding column names
methods_dict = {
    "GuaCAMOLE": "GuaCAMOLE_estimate",
    "GuaCAMOLE_eff": "GuaCAMOLE_est_eff", # Retained for potential comparison if needed
    "Bracken": "Bracken_estimate",
    "MetaPhlAn4": "mpa", # Placeholder, handled differently
}

print(f"Starting analysis for {len(parameters)} parameter sets...")

for param in parameters:

    # Check if necessary files exist before proceeding
    guac_bracken_outfile = f'results/{param}.out'
    mpa_path = f'results/{param}.mpa'
    param_file = f'parameters/{param}.csv'

    if not os.path.exists(guac_bracken_outfile):
        warnings.warn(f"Skipping {param}: Missing GuaCAMOLE/Bracken output file {guac_bracken_outfile}")
        continue
    if not os.path.exists(mpa_path) and "MetaPhlAn4" in methods_dict: # Check if MetaPhlAn4 file is needed
         warnings.warn(f"Skipping {param}: Missing MetaPhlAn4 output file {mpa_path}")
         continue
    if not os.path.exists(param_file):
        warnings.warn(f"Skipping {param}: Missing parameter file {param_file}")
        continue

    # -----------------------
    # Load GuaCAMOLE/Bracken results FIRST to check for NaNs
    # -----------------------
    try:
        tmp_df = pd.read_csv(guac_bracken_outfile, sep='\t', index_col=False)
        # Check for NaN specifically in the main GuaCAMOLE estimate column
        if tmp_df["GuaCAMOLE_estimate"].isnull().all():
            print(f"Skipping replicate {param}: GuaCAMOLE results are all NaN values.")
            skipped_replicates.append(param)
            continue # Skip to the next parameter set
    except Exception as e:
        warnings.warn(f"Skipping {param}: Error loading or checking file {guac_bracken_outfile}. Error: {e}")
        continue

    # -----------------------
    # Proceed with analysis only if GuaCAMOLE results are valid
    # -----------------------

    # 1) Load MetaPhlAn4 results (if needed and exists)
    mpa_abundances = {}
    if "MetaPhlAn4" in methods_dict:
        try:
            mpa_df = pd.read_csv(mpa_path, sep='\t', comment='#', header=None)
            # Keep only lines with s__ but not t__ => species lines
            mpa_df = mpa_df.loc[mpa_df[0].str.contains('s__') & ~mpa_df[0].str.contains('t__')]
            # Filter out rows lacking exactly 7 fields in the second column (taxid string)
            no_ids = [len([p for p in x if p]) for x in mpa_df[1].str.split('|')]
            contains_species_id = [x == 7 for x in no_ids]
            mpa_df = mpa_df.loc[contains_species_id, :]
            mpa_txids = np.array([int(i[6]) for i in mpa_df[1].str.split('|')])
            # Store {taxid -> abundance%}
            mpa_abundances = dict(zip(mpa_txids, np.array(mpa_df[2])))
        except Exception as e:
             warnings.warn(f"Could not load or parse MetaPhlAn4 results for {param}. Error: {e}. Skipping MetaPhlAn4 for this replicate.")
             # Optionally remove MetaPhlAn4 from methods_dict for this iteration if loading fails
             current_methods = {k:v for k,v in methods_dict.items() if k != "MetaPhlAn4"}
        else:
            current_methods = methods_dict.copy() # Use all methods if loaded successfully
    else:
         current_methods = {k:v for k,v in methods_dict.items() if k != "MetaPhlAn4"} # Exclude if not requested


    # -----------------------
    # 2) Set index for GuaCAMOLE/Bracken df
    # -----------------------
    # tmp_df was loaded earlier
    df = tmp_df.set_index('taxonomy_id', drop=True)

    # -----------------------
    # 3) Load true parameters
    # -----------------------
    try:
        param_df = pd.read_csv(param_file)
        # Convert #assembly_accession -> taxid
        param_df['taxid'] = [
            assemblies.loc[assembly_id, 'species_taxid']
            for assembly_id in param_df['#assembly_accession']
        ]
        param_df['gc_content'] = [
            assemblies.loc[assembly_id, 'gc_percent']
            for assembly_id in param_df['#assembly_accession']
        ]
        if param.split('_')[1] == 'medium':
            gc_medium.extend(param_df['gc_content'])
        elif param.split('_')[1] == 'extreme':
            gc_extreme.extend(param_df['gc_content'])
        elif param.split('_')[1] == 'uniform':
            gc_uniform.extend(param_df['gc_content'])
        # Dictionary of {taxid -> true_abundance}
        true_ab = dict(zip(param_df['taxid'], param_df['abundance']))
        # Sorted list of all taxids actually in this sample
        true_taxids_sorted = np.sort(list(true_ab.keys()))
    except Exception as e:
        warnings.warn(f"Skipping {param}: Error loading or processing parameter file {param_file}. Error: {e}")
        continue

    # -----------------------
    # 4) For each method, compute errors, FP, FN
    # -----------------------
    # Use current_methods which might exclude MetaPhlAn4 if it failed to load
    for method_name, column_name in current_methods.items():
        # Skip MetaPhlAn4 if it wasn't loaded successfully earlier
        if method_name == "MetaPhlAn4" and not mpa_abundances:
            continue

        # ab_true = vector of true abundances for each taxid in sorted order
        ab_true_vec = np.array([true_ab[t] for t in true_taxids_sorted])

        if method_name == "MetaPhlAn4":
            # Convert from percentage to fraction
            ab_method_vec = np.array([
                mpa_abundances[t] if t in mpa_abundances else 0
                for t in true_taxids_sorted
            ]) / 100.0
        else:
            # e.g. GuaCAMOLE, Bracken, etc.
            # Ensure index access is safe (check if taxid exists)
            ab_method_vec = np.array([
                df.loc[t, column_name] if t in df.index and column_name in df.columns else 0
                for t in true_taxids_sorted
            ])

        # Relative error calculation with check for zero true abundance
        # re_method = (ab_method_vec - ab_true_vec) / ab_true_vec # Old way
        re_method = np.full_like(ab_true_vec, np.nan, dtype=np.float64) # Initialize with NaN
        non_zero_true_mask = ab_true_vec != 0
        re_method[non_zero_true_mask] = (
            ab_method_vec[non_zero_true_mask] - ab_true_vec[non_zero_true_mask]
            ) / ab_true_vec[non_zero_true_mask]


        # Store entry-by-entry detail (only if re_method is not all NaN)
        if not np.isnan(re_method).all():
            txid_error.extend(re_method) # Contains NaNs where true abundance was 0
            abundance_list.extend(ab_method_vec)
            method_txid.extend([method_name] * len(re_method))
            sample_txid.extend([param[:-2]] * len(re_method)) # Store sample group identifier
            replicate_txid.extend([param] * len(re_method)) # Store full replicate id
            txid_list.extend(true_taxids_sorted)


        # Summarize error for the sample (Mean Absolute Relative Error on detected & true taxa)
        # Mask for taxa that are truly present AND detected by the method
        valid_comparison_mask = (ab_true_vec > 0) & (ab_method_vec > 0)

        re_method_filtered = re_method[valid_comparison_mask]

        if len(re_method_filtered) > 0:
            # Use nanmean just in case (though NaNs from true=0 were excluded by mask)
            mean_abs_err = np.nanmean(np.abs(re_method_filtered))
        else:
            # Handle case where no common taxa are detected > 0
            mean_abs_err = np.nan

        error.append(mean_abs_err)
        methods.append(method_name)
        samples_list.append(param[:-2]) # Store sample group identifier
        replicates.append(param) # Store full replicate id

        # (FP/FN calculation remains the same as your original code)
        # 1) taxids truly in sample with ab>0
        truth_positive_set = set(
            t for t, val in true_ab.items()
            if val > 0
        )

        if method_name == "MetaPhlAn4":
            method_positive_set = set(
                t for t, val in mpa_abundances.items()
                if val > 0
            )
        else:
            method_positive_set = set(
                t for t in df.index
                if t in df[column_name] and pd.notna(df.loc[t, column_name]) and df.loc[t, column_name] > 0
            ) # Added notna check

        fn_count = len(truth_positive_set - method_positive_set)
        all_taxids_in_sample = set(true_ab.keys())
        fp_count = len(method_positive_set - truth_positive_set) # FP defined as detected > 0 when true = 0

        f_positives.append(fp_count)
        f_negatives.append(fn_count)


# --- Report Skipped Replicates ---
print("\n--------------------------------------")
print(f"Analysis complete.")
print(f"Total replicates skipped due to GuaCAMOLE NaN values: {len(skipped_replicates)}")
if skipped_replicates:
    print("Skipped replicate identifiers:")
    for rep_id in skipped_replicates:
        print(f"  - {rep_id}")
print("--------------------------------------\n")
# ----------------------------------


# -----------------------
# Build err_df with FP, FN
# -----------------------
err_df = pd.DataFrame({
    'error': [x*100 for x in error],
    'method': methods,
    'sample_raw': samples_list, # Contains raw group identifier like 'param_uniform_400'
    'replicate': replicates,    # Contains full replicate identifier like 'param_uniform_400_rep_1'
    'false_positive': f_positives,
    'false_negative': f_negatives,
})

# Derive GC_dist & n_species from the raw sample group identifier
# Assuming format like 'param_dist_norg' stored in samples_list -> sample_raw
err_df['GC_dist'] = err_df['sample_raw'].apply(lambda x: x.split('_')[1])
err_df['n_species'] = err_df['sample_raw'].apply(lambda x: x.split('_')[2])
# Create the clean sample group identifier (e.g., 'uniform_400') for plotting
err_df['sample_group'] = err_df['sample_raw'].apply(lambda x: '_'.join(x.split('_')[1:3]))


# And if needed, the per-taxon info frame
txid_df = pd.DataFrame({
    'abundance': abundance_list,
    'method': method_txid,
    # 'sample' column here should also use the group identifier if needed for consistency
    'sample_group_raw': sample_txid, # Contains raw group identifier like 'param_uniform_400'
    'txid': txid_list,
    'error': txid_error, # Contains NaNs where true=0
    'replicate': replicate_txid # Contains full replicate id
})
# Add the clean sample group identifier to txid_df as well
txid_df['sample_group'] = txid_df['sample_group_raw'].apply(lambda x: '_'.join(x.split('_')[1:3]))

# %%
#plot_order = ["GuaCAMOLE", "GuaCAMOLE_eff", "Bracken", "MetaPhlAn4", "GuaCAMOLE_eff"]
#final_legend_order = ["GuaCAMOLE", , "Bracken", "MetaPhlAn4"]
plot_order = ["GuaCAMOLE", "Bracken", "MetaPhlAn4"]
final_legend_order = ["GuaCAMOLE", "Bracken", "MetaPhlAn4"]

# Example of recreating x_order for plots using the correct 'sample_group' column
# (Keep your existing logic for creating x_order)
unique_sample_groups = sorted(err_df['sample_group'].unique())
gc_dists_plot = sorted(list(set(err_df['GC_dist'])))
n_species_plot = sorted(list(set(err_df['n_species'].astype(int))))
x_order_intended = ['_'.join([gc, str(species)]) for gc in gc_dists_plot for species in n_species_plot]
x_order = [s for s in x_order_intended if s in unique_sample_groups]

print("Plotting overall error...")
plt.figure(figsize=(6, 4))
# Filter data for plot (as provided by user)
err_df_filtered_for_plot = err_df[
    (~(err_df['replicate'] == 'param_medium_5_rep_0') &
    (~(err_df['method'] == 'GuaCAMOLE_eff')))
].copy()

# Plot using specified order and palette
ax_err = sns.stripplot(x='sample_group', y='error', hue='method',
                       data=err_df_filtered_for_plot,
                       order=x_order,
                       hue_order=plot_order, # Control plot/dodge order (GuaCAMOLE last)
                       palette=method_colors, # Use defined colors
                       dodge=True)

# --- Add Vertical Separator Lines ---
num_categories = len(x_order)
# Calculate positions for lines separating the 3 groups
# Assumes x_order has categories grouped by GC strategy
# Line after the first third, line after the second third
pos1 = num_categories / 3 - 0.5
pos2 = 2 * num_categories / 3 - 0.5

# Check if positions are valid (only draw if num_categories is divisible by 3 or handled appropriately)
ax_err.axvline(x=pos1, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0) # zorder=0 sends it behind points
ax_err.axvline(x=pos2, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0)

# Set custom xticklabels
new_xticklabels_err = [label.split('_')[-1] for label in x_order]
ax_err.set_xticklabels(new_xticklabels_err)
ax_err.tick_params(axis='x', rotation=0) # Reset rotation if desired, or set to 45

# Manually rebuild legend with desired order and no frame
handles, labels = ax_err.get_legend_handles_labels()
label_handle_dict_err = dict(zip(labels, handles))
ordered_handles_err = [label_handle_dict_err[label] for label in final_legend_order if label in label_handle_dict_err]
ordered_labels_err = [label for label in final_legend_order if label in label_handle_dict_err]
ax_err.legend(ordered_handles_err, ordered_labels_err, frameon=False) # Example placement outside

plt.xlabel("Number of Taxa in Sample") # Simplified label
plt.ylabel("Mean Absolute Relative Error (%)")
plt.tight_layout() # Adjust layout to make space for legend outside
plt.savefig('err_strip_plot_wo_eff.pdf',
            transparent=True, bbox_inches='tight')

# print some stats:
g = err_df.groupby(['method'])['false_positive'].mean()

# %%
print("Plotting False Positives...")
plt.figure(figsize=(6, 4))
# Plot using specified order and palette
ax_fp = sns.stripplot(x='sample_group', y='false_positive', hue='method', data=err_df,
                      order=x_order,
                      hue_order=plot_order, # Control plot/dodge order
                      palette=method_colors, # Use defined colors
                      dodge=True)

# Set custom xticklabels (use labels generated previously)
ax_fp.axvline(x=pos1, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0) # zorder=0 sends it behind points
ax_fp.axvline(x=pos2, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0)
ax_fp.set_xticklabels(new_xticklabels_err)
ax_fp.tick_params(axis='x', rotation=0) # Reset rotation if desired, or set to 45

# Manually rebuild legend
handles_fp, labels_fp = ax_fp.get_legend_handles_labels()
label_handle_dict_fp = dict(zip(labels_fp, handles_fp))
ordered_handles_fp = [label_handle_dict_fp[label] for label in final_legend_order if label in label_handle_dict_fp]
ordered_labels_fp = [label for label in final_legend_order if label in label_handle_dict_fp]
ax_fp.legend(ordered_handles_fp, ordered_labels_fp, frameon=False)

plt.xlabel("Number of Taxa in Sample") # Simplified label
plt.ylabel("False Positive Taxa")
plt.tight_layout() # Adjust layout to make space for legend outside
plt.savefig('fp_strip_plot.pdf',
            transparent=True, bbox_inches='tight')

print("Plotting False Negatives...")
plt.figure(figsize=(6, 4))
# Plot using specified order and palette
ax_fn = sns.stripplot(x='sample_group', y='false_negative', hue='method', data=err_df,
                      order=x_order,
                      hue_order=plot_order, # Control plot/dodge order
                      palette=method_colors, # Use defined colors
                      dodge=True)

# Set custom xticklabels
ax_fn.axvline(x=pos1, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0) # zorder=0 sends it behind points
ax_fn.axvline(x=pos2, color='grey', linestyle='-', linewidth=1, alpha=0.7, zorder=0)

ax_fn.set_xticklabels(new_xticklabels_err)
ax_fn.tick_params(axis='x', rotation=0) # Reset rotation if desired, or set to 45

# Manually rebuild legend
handles_fn, labels_fn = ax_fn.get_legend_handles_labels()
label_handle_dict_fn = dict(zip(labels_fn, handles_fn))
ordered_handles_fn = [label_handle_dict_fn[label] for label in final_legend_order if label in label_handle_dict_fn]
ordered_labels_fn = [label for label in final_legend_order if label in label_handle_dict_fn]
ax_fn.legend(ordered_handles_fn, ordered_labels_fn, frameon=False)

plt.xlabel("Number of Taxa in Sample") # Simplified label
plt.ylabel("False Negative Taxa")
plt.tight_layout() # Adjust layout to make space for legend outside
plt.savefig('fn_strip_plot.pdf',
            transparent=True, bbox_inches='tight')

