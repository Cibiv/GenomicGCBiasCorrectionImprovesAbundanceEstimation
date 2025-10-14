#!/usr/bin/env python3
"""
bulk_iss_generate_selective.py

Usage:
  python bulk_iss_generate_selective.py

Description:
  - Iterates over all CSV parameter files in the 'parameters/' folder.
  - Each has columns: '#assembly_accession' and 'abundance'.
  - For each parameter file:
    1. Collect the needed assemblies.
    2. Locate each assembly's .fna file in 'assemblies/'.
    3. Create a small abundance file with no header for ISS.
    4. Call ISS with exactly those files in '-g'.

Example param file row:
  #assembly_accession,abundances
  GCF_000368265.1,0.2
  GCF_000368266.1,0.8
"""

import os
import subprocess
import pandas as pd

def main():
    param_dir = "parameters"
    assemblies_dir = "assemblies_merged"

    # Ensure directories exist
    if not os.path.isdir(param_dir):
        print(f"Error: '{param_dir}' folder not found.")
        return
    if not os.path.isdir(assemblies_dir):
        print(f"Error: '{assemblies_dir}' folder not found.")
        return

    # List all CSV parameter files in 'parameters/'.
    param_files = [
        f for f in os.listdir(param_dir)
        if os.path.isfile(os.path.join(param_dir, f))
           and f.endswith('.csv')
           and f.startswith('param_')
    ]

    # We'll also make a set of all .fna files we have in 'assemblies' to confirm existence
    available_fnas = set(os.listdir(assemblies_dir))

    for param_file in param_files:
        param_path = os.path.join(param_dir, param_file)

        # Derive a prefix for naming the output from the param file name
        # e.g., param_medium_100_rep_0.csv -> param_medium_100_rep_0
        prefix, _ = os.path.splitext(param_file)

        distribution = param_file.split('_')[1]
        if distribution == "extreme":
            nreads = "45M"
        else:
            nreads = "40M"

        # Read the parameter CSV
        try:
            df = pd.read_csv(param_path, sep=None, engine='python')  # auto-detect delimiter
        except Exception as e:
            print(f"Skipping {param_file}, could not parse as CSV: {e}")
            continue

        # Check the required columns
        if '#assembly_accession' not in df.columns or 'abundance' not in df.columns:
            print(f"Skipping {param_file}, missing required columns.")
            continue

        # Build the list of needed assemblies
        # e.g. GCF_000368265.1 from #assembly_accession
        needed_assemblies = df['#assembly_accession'].tolist()
        needed_abundances = df['abundance'].tolist()

        # Double-check we have .fna files for them
        fna_filepaths = []
        for asm in needed_assemblies:
            # Typically the .fna might be named GCF_000368265.1.fna
            # If that's the pattern, we do:
            candidate = asm + ".fna"
            if candidate in available_fnas:
                fna_filepaths.append(os.path.join(assemblies_dir, candidate))
            else:
                # If there's no direct match, we might look for partial matches, or skip
                print(f"WARNING: no .fna file found for assembly {asm}")
                # You could skip or handle differently
                # For now let's just skip
                continue

        # If we have no .fna files, skip
        if not fna_filepaths:
            print(f"No matching assemblies found for {param_file}, skipping.")
            continue

        # Next, create a temp abundance file (two columns, no header)
        abundance_temp = f"{prefix}_abundance.temp.txt"
        with open(abundance_temp, 'w') as temp_f:
            for asm, ab in zip(needed_assemblies, needed_abundances):
                temp_f.write(f"{asm}\t{ab}\n")

        # Build the iss command
        # Example:
        #   iss generate -g <file1.fna> <file2.fna> ... \
        #       --model hiseq --cpus 10 \
        #       --abundance_file <temp_file> \
        #       -o <prefix>
        if not os.path.exists(os.path.join('libraries/', prefix)):
            os.mkdir(os.path.join('libraries/', prefix))

        out_prefix = os.path.join('libraries/', prefix, prefix)
        cmd = [
            "iss", "generate",
            "--model", "hiseq",
            "--cpus", "80",
            "--n_reads", nreads,
            "--gc_bias",
            "--abundance_file", abundance_temp,
            "-o", out_prefix,
            "-g",  # after this we list all the .fna paths
        ] + fna_filepaths

        print(f"\n=== Running ISS for: {param_file} ===")
        print("Command:", " ".join(cmd))

        try:
            subprocess.run(cmd, check=True)
            print(f"Done generating reads for {param_file}. Output prefix: {prefix}")
        except subprocess.CalledProcessError as e:
            print(f"Error while running ISS for {param_file}: {e}")
        finally:
            # Clean up the abundance file
            if os.path.exists(abundance_temp):
                os.remove(abundance_temp)
            subprocess.run('rm ' + out_prefix + '*.vcf', shell=True)

        r1 = f"{out_prefix}_R1.fastq"
        r2 = f"{out_prefix}_R2.fastq"

        # 3) Spawn gzip in the background
        # We can pass multiple files to gzip at once
        if os.path.exists(r1) and os.path.exists(r2):
            gzip_cmd = ["gzip", r1, r2]
            print(f"Compressing FASTQs in the background: {' '.join(gzip_cmd)}")
            subprocess.Popen(gzip_cmd)  # does NOT block
        else:
            print("WARNING: One or both FASTQ files are missing, skipping compression.")

if __name__ == "__main__":
    main()

