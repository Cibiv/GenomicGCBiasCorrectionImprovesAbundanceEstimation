# Description

This is the companion repository for our publication "Genomic GC bias correction improves species abundance estimation from metagenomic data" where discuss and evaluate our GC-bias correction algorithm [GuaCAMOLE](https://github.com/Cibiv/GuaCAMOLE) using a large number of datasets. This repository contains all necessary raw data and scripts to reproduce the main analyses and figures shown in our publication.

For the [GuaCAMOLE](https://github.com/Cibiv/GuaCAMOLE) algorithm, please go to [https://github.com/Cibiv/GuaCAMOLE](https://github.com/Cibiv/GuaCAMOLE).

# Contents

## Fig. 2B, simulation study (``simulation/``)

The script ``simulate_communities.py`` simulates communities with different number of taxa and genomic GC composition as described in the manuscript. Each simulated community is represented by a file in ``parameters/`` which contains the RefSeq identifier of the genome assembly used to represent each taxon and the taxon's abundance in the community. All assemblies are downloaded into ``assemblies_merged`` in a format suitable for read simulation. The script ``simulate_reads.py`` uses a modified version of ISS to simulate sequencing reads for each libraries and places the simulated libraries in `libraries/`. `evaluate.py` reads the output for GuaCAMOLE and MetaPhlAn4 from `results/` and reproduces the simulation results shown in Fig. 2. Due to size constraints, this repository contains the community descriptions and GuaCAMOLE/MetaPhlAn4 results, but not the genome assemblies (these are available from RefSeq) and simulated reads.

## Fig. 3, data from Tourlousse et al., 2024 (``tourlousse/``)

The file `sra_ids.txt` lists the SRA ids of the samples of Tourlousse et al. (2021) used to evaluate GuaCAMOLE. The folders `guacamole_results`, `motus_results`, `singlem_results` and `sylph_results` contain the output produced by these tools for the data of Tourlousse et al. `evaluate_efficiencies.py` reads the GuaCAMOLE results and reproduces Fig. 3A. `evaluate.py` the output of all tools and reproduces panels Fig. 3 B-E.

## Fig. 4, colorectal cancer datasets (``crc/``)

The files `gupta2020_s3.csv`, `gupta2020_s6.csv`, and `murovec2024_s2.csv`, `PRJDB4176.csv`, `PRJDB4176.attributes.tab` and `PRJDB4176.sample_disease_status.csv` contain the supplemental tables of Gupta et al. (2020) and Murovec at al. (2024) from which the list of samples (`samples.yaml`) used was derived as described in the publication. The RStudio Notebook `evaluate.Rmd` reads these samples lista and the raw GuaCAMOLE results from `guacamole_results` and perform the clustering and analysis as described in the publication. It produces the plots shwon in Fig. 4 and additional generates listing all individual abundance estimates (`abundances.csv.bz2`), and all studies including their cluster assignment (`studies.csv.bz2`).

## Fig. S1, data from Mori et al., 2023 (``mori/``)

The file `samples.csv` list the samples from Mori et al. (2023) comparing abundance estimates across sequencing platforms and DNA extraction methods. `evaluate.Rmd` reads the raw GuaCAMOLE results for these samples and reproduces the panels shown in Fig. S1.
