#!/usr/bin/env python3

import os
import numpy as np
import pandas as pd
import itertools
import urllib.request
import random
import gzip
import shutil
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

###########################
# Utility: Merge contigs into one record
###########################
def merge_contigs_and_rename(in_fna_path, out_fna_path, assembly_id):
    """
    Reads all FASTA records from in_fna_path (which may be .gz or uncompressed),
    merges them into one single sequence,
    and writes out to out_fna_path with a single record named assembly_id.
    """
    merged_seq_str = ""
    # Open .gz or uncompressed
    if in_fna_path.endswith(".gz"):
        handle = gzip.open(in_fna_path, "rt")  # text mode
    else:
        handle = open(in_fna_path, "r")

    with handle:
        records = list(SeqIO.parse(handle, "fasta"))
        if not records:
            raise ValueError(f"No valid FASTA records in {in_fna_path}")
        # Concatenate all contigs
        merged_seq_str = "".join(str(rec.seq) for rec in records)

    # Write single-record FASTA with the ID = assembly_id
    merged_record = SeqRecord(
        Seq(merged_seq_str),
        id=assembly_id,
        description=""
    )
    with open(out_fna_path, "w") as out_f:
        SeqIO.write(merged_record, out_f, "fasta")

###########################
# Utility: Ensure local assembly
###########################
def ensure_assembly_exists(assembly_id, ftp_path, assemblies_dir):
    """
    Checks if assemblies_dir/assembly_id.fna already exists as a single-record FASTA.
    If not, attempts to download from ftp_path, merges contigs, renames the record.
    If success, returns True; if fails, returns False.
    ftp_path is expected to be something like:
        'ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/007/745.1'
    from which we build the file name: assembly_id + '_genomic.fna.gz'
    and then handle merging.
    """
    print('calling ensure assembly exists fucntion')
    out_fna = os.path.join(assemblies_dir, f"{assembly_id}.fna")
    if os.path.isfile(out_fna):
        # Already have it
        return True

    # Attempt download
    assembly_name = ftp_path.split('/')[-1] + '_genomic.fna.gz'
    download_url = ftp_path + "/" + assembly_name
    file_name_zip = assembly_id + '.fna.gz'
    file_name = assembly_id + '.fna'

    # If already partially downloaded, skip re-download
    if not os.path.exists(out_fna):
        try:
            print(f"Downloading {download_url} -> {file_name_zip}")
            urllib.request.urlretrieve(download_url, os.path.join(assemblies_dir, file_name_zip))
        except Exception as e:
            print(f"Failed to download {download_url}: {e}")
            return False

    # Now merge contigs and unzip
    try:
        merge_contigs_and_rename(os.path.join(assemblies_dir, file_name_zip), out_fna, assembly_id)
    except Exception as e:
        print(f"Failed to merge contigs for {assembly_id}: {e}")
        return False

    # If successful, optionally remove the .gz file to save space:
    # os.remove(tmp_download)

    return True


###########################
# Distributions
###########################
def make_medium_distribution(nbins=10, center=None, sigma=None):
    if center is None:
        center = (nbins - 1) / 2.0
    if sigma is None:
        sigma = nbins / 6.0
    x = np.arange(nbins)
    bell = np.exp(-0.5 * ((x - center) / sigma) ** 2)
    bell /= bell.sum()
    return bell

def make_extreme_distribution(nbins=10, center=None, sigma=None):
    bell = make_medium_distribution(nbins, center, sigma)
    inverted = 1.0 - bell / max(bell)
    inverted[inverted < 0] = 0
    inverted /= inverted.sum()
    return inverted

def sample_binned_gc(
        df,
        num_samples=50,
        gc_col='gc_percent',
        bin_edges=None,
        distribution='uniform',
        random_state=None
):
    if bin_edges is None:
        bin_edges = [25,30,35,40,45,50,55,60,65,70,75]

    df = df.dropna(subset=[gc_col]).copy()
    df['gc_bin'] = pd.cut(df[gc_col], bins=bin_edges, include_lowest=True)

    bin_labels = df['gc_bin'].cat.categories
    num_bins = len(bin_labels)

    if distribution == 'uniform':
        dist = np.ones(num_bins) / num_bins
    elif distribution == 'medium':
        dist = make_medium_distribution(num_bins)
    elif distribution == 'extreme':
        dist = make_extreme_distribution(num_bins)
    else:
        raise ValueError("distribution must be 'uniform', 'medium', or 'extreme'")

    if not np.isclose(dist.sum(), 1.0):
        dist /= dist.sum()

    bin_samples = (dist * num_samples).round().astype(int)
    diff = num_samples - bin_samples.sum()
    while diff != 0:
        if diff > 0:
            idx = np.random.choice(num_bins)
            bin_samples[idx] += 1
            diff -= 1
        else:
            idx_candidates = np.where(bin_samples > 0)[0]
            if len(idx_candidates) == 0:
                break
            idx = np.random.choice(idx_candidates)
            bin_samples[idx] -= 1
            diff += 1

    grouped = df.groupby('gc_bin')
    sampled_rows = []
    for bin_label, n_samp in zip(bin_labels, bin_samples):
        if n_samp <= 0:
            continue
        group_df = grouped.get_group(bin_label)

        # We enforce distinct species in that bin
        grouped_by_species = group_df.groupby('species_taxid')
        unique_species = list(grouped_by_species.groups.keys())
        n_pick = min(n_samp, len(unique_species))

        chosen_species = np.random.choice(unique_species, size=n_pick, replace=False)
        chosen_rows_bin = []
        for s in chosen_species:
            species_df = grouped_by_species.get_group(s)
            chosen_one = species_df.sample(n=1, replace=False, random_state=random_state)
            chosen_rows_bin.append(chosen_one)
        chosen = pd.concat(chosen_rows_bin, ignore_index=True)
        sampled_rows.append(chosen)

    final_sample = pd.concat(sampled_rows).reset_index(drop=True)
    return final_sample


###########################
# MAIN SCRIPT
###########################
def main():

    # Load your DataFrame
    df = pd.read_csv(
        'assembly_summary.txt.bz2',
        sep='\t', 
        header=0
    )

    # We'll store parameter CSVs here
    out_dir = "parameters"
    os.makedirs(out_dir, exist_ok=True)

    assemblies_dir = "assemblies_merged"  # Where we keep single-record .fna
    os.makedirs(assemblies_dir, exist_ok=True)

    organism_numbers = [5, 10, 50, 100, 400]
    replications = 5
    distributions = ['uniform', 'medium', 'extreme']

    for dist, n_org, rep in itertools.product(distributions, organism_numbers, range(replications)):
        print(f"\n=== Sampling distribution={dist}, n_org={n_org}, replicate={rep} ===")

        # We'll do a small loop that tries re-sampling if an assembly can't be acquired
        max_tries = 50  # avoid infinite loop
        for attempt in range(max_tries):
            sampled_df = sample_binned_gc(df, num_samples=n_org, distribution=dist)
            # Check each row's assembly
            all_good = True
            for idx, row in sampled_df.iterrows():
                asm_acc = row['#assembly_accession']
                ftp_path = row['ftp_path']
                # Ensure we have the .fna in assemblies_merged
                if not ensure_assembly_exists(asm_acc, ftp_path, assemblies_dir):
                    all_good = False
                    print(f"Assembly {asm_acc} missing or download failed. Retrying sampling.")
                    break

            if all_good:
                # All assemblies exist, so we finalize
                abundances = np.random.lognormal(size=sampled_df.shape[0])
                sampled_df['abundance'] = abundances / abundances.sum()

                out_path = os.path.join(out_dir, f'param_{dist}_{n_org}_rep_{rep}.csv')
                # The user wants columns: #assembly_accession, abundance
                # plus anything else you want. We'll do species_taxid, etc. if needed
                final_df = sampled_df[['#assembly_accession', 'abundance']].copy()
                final_df.to_csv(out_path, index=False)
                print(f"Wrote {out_path} with {len(final_df)} assemblies.")
                break
        else:
            print(f"Failed to get a valid sample after {max_tries} tries. Skipping {dist}, {n_org}, {rep}.")

if __name__ == "__main__":
    main()

