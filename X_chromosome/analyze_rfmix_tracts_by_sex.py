#!/usr/bin/env python3
"""
Analyze RFMix output on X chromosome and calculate ancestry tract lengths
separated by sex and by sex x ancestry.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import gzip
import os
import argparse

 
# I/O
 

def read_hmm_file(hmm_file_path):
    print(f"Reading HMM file: {hmm_file_path}")
    hmm = pd.read_csv(hmm_file_path, sep=r'\s+', header=None)
    hmm.rename(columns={0: 'chr', 1: 'pos', 6: 'morgan'}, inplace=True)
    return hmm


def read_rfmix_file(path):
    if path.endswith(".gz"):
        with gzip.open(path, "rt") as f:
            data = [list(map(int, line.split())) for line in f if line.strip()]
    else:
        with open(path) as f:
            data = [list(map(int, line.split())) for line in f if line.strip()]
    return np.array(data)


def read_sex_ploidy(path):
    """
    Returns a list mapping haplotype index -> sex
    """
    df = pd.read_csv(path, sep=r'\s+', header=None)
    df.columns = ["sample", "sex", "aux"]

    hap_sex = []
    for _, row in df.iterrows():
        if row.sex == 2:      # female
            hap_sex.extend(["female", "female"])
        elif row.sex == 1:    # male
            hap_sex.append("male")
        else:
            raise ValueError(f"Unknown sex code: {row.sex}")

    return hap_sex


 
# Tract calculation
 

def calculate_tracts(rfmix, morgan, hap_sex):
    """
    Returns:
      tracts_by_sex
      tracts_by_sex_and_ancestry
    """

    n_snps, n_haps = rfmix.shape

    morgan = morgan[:n_snps]
    cum_morgan = np.cumsum([0] + list(morgan[1:]))

    tracts_by_sex = {"male": [], "female": []}
    tracts_by_sex_anc = {
        ("male", 1): [], ("male", 2): [],
        ("female", 1): [], ("female", 2): []
    }

    for h in range(n_haps):
        ancestry = rfmix[:, h]
        sex = hap_sex[h]

        start = 0
        current = ancestry[0]

        for i in range(1, n_snps):
            if ancestry[i] != current:
                length_cm = (cum_morgan[i-1] - cum_morgan[start]) * 100
                tracts_by_sex[sex].append(length_cm)
                tracts_by_sex_anc[(sex, current)].append(length_cm)

                start = i
                current = ancestry[i]

        # last tract
        length_cm = (cum_morgan[n_snps-1] - cum_morgan[start]) * 100
        tracts_by_sex[sex].append(length_cm)
        tracts_by_sex_anc[(sex, current)].append(length_cm)

    return tracts_by_sex, tracts_by_sex_anc


 
# Plotting
 

def plot_sex(tracts, out):
    plt.figure(figsize=(7,6))
    bins = np.logspace(0, 3, 50)

    for sex, color in [("male", "blue"), ("female", "red")]:
        counts, edges = np.histogram(tracts[sex], bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2
        plt.scatter(centers[counts>0], counts[counts>0], label=sex, color=color)

    plt.yscale("log")
    plt.xlabel("Tract length (cM)")
    plt.ylabel("Number of tracts")
    plt.title("X ancestry tract lengths by sex")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


def plot_sex_ancestry(tracts, out):
    plt.figure(figsize=(7,6))
    bins = np.logspace(0, 3, 50)

    styles = {
        ("male", 1): ("blue", "Male A1"),
        ("male", 2): ("cyan", "Male A2"),
        ("female", 1): ("red", "Female A1"),
        ("female", 2): ("orange", "Female A2"),
    }

    for key, (color, label) in styles.items():
        counts, edges = np.histogram(tracts[key], bins=bins)
        centers = (edges[:-1] + edges[1:]) / 2
        plt.scatter(centers[counts>0], counts[counts>0], label=label, color=color)

    plt.yscale("log")
    plt.xlabel("Tract length (cM)")
    plt.ylabel("Number of tracts")
    plt.title("X ancestry tract lengths by sex and ancestry")
    plt.legend()
    plt.tight_layout()
    plt.savefig(out, dpi=300)
    plt.close()


 
# Main
 

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hmm-file", required=True)
    ap.add_argument("--rfmix-pattern", required=True)
    ap.add_argument("--sex-ploidy", required=True)
    ap.add_argument("--output-prefix", default="rfmix")
    ap.add_argument("--chr-start", type=int, default=23)
    ap.add_argument("--chr-end", type=int, default=23)
    args = ap.parse_args()

    hmm = read_hmm_file(args.hmm_file)
    hap_sex = read_sex_ploidy(args.sex_ploidy)

    tracts_sex_all = {"male": [], "female": []}
    tracts_sex_anc_all = {
        ("male",1):[],("male",2):[],
        ("female",1):[],("female",2):[]
    }

    for chrn in range(args.chr_start, args.chr_end + 1):
        rfmix_file = args.rfmix_pattern.replace("${chr}", str(chrn))
        if not os.path.exists(rfmix_file):
            continue

        rfmix = read_rfmix_file(rfmix_file)
        morgan = hmm[hmm.chr == chrn]["morgan"].values

        t_sex, t_sa = calculate_tracts(rfmix, morgan, hap_sex)

        for k in tracts_sex_all:
            tracts_sex_all[k].extend(t_sex[k])
        for k in tracts_sex_anc_all:
            tracts_sex_anc_all[k].extend(t_sa[k])

    # write outputs
    with open(f"{args.output_prefix}_tracts_by_sex.txt","w") as f:
        f.write("sex\ttract_length_cm\n")
        for sex in tracts_sex_all:
            for l in tracts_sex_all[sex]:
                f.write(f"{sex}\t{l:.6f}\n")

    with open(f"{args.output_prefix}_tracts_by_sex_and_ancestry.txt","w") as f:
        f.write("sex\tancestry\ttract_length_cm\n")
        for (sex, anc), vals in tracts_sex_anc_all.items():
            for l in vals:
                f.write(f"{sex}\t{anc}\t{l:.6f}\n")

    plot_sex(tracts_sex_all, f"{args.output_prefix}_by_sex.png")
    plot_sex_ancestry(tracts_sex_anc_all, f"{args.output_prefix}_by_sex_and_ancestry.png")


if __name__ == "__main__":
    main()
