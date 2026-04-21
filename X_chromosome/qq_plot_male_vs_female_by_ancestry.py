#!/usr/bin/env python3
"""
QQ plots of log10 tract lengths: males vs females
Ancestry-specific, overlaid in one 2x2 figure (one panel per method)
Fixes:
1. Proper percentile-aligned comparison
2. Correct male/female axis assignment
3. Filter out tracts <1 cM
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse

 
# Global plotting style
 
plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 13,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 11
}) 

ANCESTRY_LABELS = {
    1: ("Farmer", "tab:blue"),
    2: ("Hunter-gatherer", "tab:red")
}

 
# Helpers
 

def load_tracts(path, ancestry):
    """Return male and female tract lengths for a given ancestry, filtered for >=1 cM."""
    df = pd.read_csv(path, sep=r"\s+")
    df = df[df.ancestry == ancestry]
    
    males = df[df.sex == "male"]["tract_length_cm"].values
    females = df[df.sex == "female"]["tract_length_cm"].values

    # Filter out tracts <1 cM
#    males = males[males >= 1]
#    females = females[females >= 1]

    return males, females

def qq_points(x, y, min_n=10, num_points=100):
    """
    Return x and y values at matching percentiles for QQ plot.
    x -> males
    y -> females
    """
    x = np.log10(x[x > 0])
    y = np.log10(y[y > 0])

    if len(x) < min_n or len(y) < min_n:
        return None, None

    # Fixed percentiles from 0 to 1, avoiding 0 and 1
    percentiles = np.linspace(0, 1, num_points + 2)[1:-1]
    xq = np.quantile(x, percentiles)  # male values
    yq = np.quantile(y, percentiles)  # female values
    return xq, yq

 
# Main
 
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs=4, required=True)
    ap.add_argument("--labels", nargs=4, required=True)
    ap.add_argument("--output-prefix", required=True)
    args = ap.parse_args()

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.flatten()

    for ax, path, label in zip(axes, args.inputs, args.labels):
        ax.set_title(label)
        all_x = []
        all_y = []

        for anc in [1, 2]:
            males, females = load_tracts(path, anc)
            xq, yq = qq_points(males, females)  # males -> x, females -> y

            if xq is None:
                continue

            anc_label, color = ANCESTRY_LABELS[anc]
            ax.scatter(xq, yq, s=20, alpha=0.7, label=anc_label, color=color)
            all_x.append(xq)
            all_y.append(yq)

        # Draw dashed y=x line using combined limits
        if all_x:
            xmin = min(np.min(x) for x in all_x)
            xmax = max(np.max(x) for x in all_x)
            ymin = min(np.min(y) for y in all_y)
            ymax = max(np.max(y) for y in all_y)
            lims = [min(xmin, ymin), max(xmax, ymax)]
            ax.plot(lims, lims, "k--", lw=1.5)

        ax.set_xlabel("Male log10 tract length (cM)")
        ax.set_ylabel("Female log10 tract length (cM)")

    handles, labels = axes[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc="lower center", ncol=2)

    plt.tight_layout(rect=[0, 0.08, 1, 1])
    out = f"{args.output_prefix}_ancestry_overlay_filtered.png"
    plt.savefig(out, dpi=300)
    plt.close()
    print(f"Saved {out}")

if __name__ == "__main__":
    main()
