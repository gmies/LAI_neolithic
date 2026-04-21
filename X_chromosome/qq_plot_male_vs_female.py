#!/usr/bin/env python3
"""
QQ plots of log10 tract lengths: males (x) vs females (y)
One 2x2 panel (one per method)
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

 
# Helpers
 

def load_tracts(path):
    df = pd.read_csv(path, sep=r"\s+")
    males = df[df.sex == "male"]["tract_length_cm"].values
    females = df[df.sex == "female"]["tract_length_cm"].values
    return males, females


def qq_points(x, y, min_n=10):
    x = np.log10(x[x > 0])
    y = np.log10(y[y > 0])

    n = min(len(x), len(y))
    if n < min_n:
        return None, None

    q = np.linspace(0, 1, n + 1)[1:-1]
    return np.quantile(x, q), np.quantile(y, q)

 
# Main
 

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inputs", nargs=4, required=True)
    ap.add_argument("--labels", nargs=4, required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    fig, axes = plt.subplots(2, 2, figsize=(10, 10))
    axes = axes.flatten()

    for ax, path, label in zip(axes, args.inputs, args.labels):
        males, females = load_tracts(path)
        xq, yq = qq_points(males, females)

        ax.set_title(label)

        if xq is None:
            ax.text(0.5, 0.5, "Insufficient data",
                    ha="center", va="center", transform=ax.transAxes)
            ax.set_xlabel("Male log10 tract length (cM)")
            ax.set_ylabel("Female log10 tract length (cM)")
            continue

        ax.scatter(xq, yq, s=20, alpha=0.7)

        lims = [
            min(xq.min(), yq.min()),
            max(xq.max(), yq.max())
        ]
        ax.plot(lims, lims, "k--", lw=1.5)

        ax.set_xlabel("Male log10 tract length (cM)")
        ax.set_ylabel("Female log10 tract length (cM)")

    plt.tight_layout()
    plt.savefig(args.output, dpi=300)
    plt.close()

    print(f"Saved {args.output}")


if __name__ == "__main__":
    main()
