#!/usr/bin/env python3

import numpy as np
import matplotlib
matplotlib.use('Agg')  # Add this line to prevent GUI backend issues
import matplotlib.pyplot as plt
import sys
import os


def detect_elbow(bit_scores):
    """
    Use the geometric (knee) method to detect the elbow point in a sorted list of scores.
    """
    scores = np.array(sorted(bit_scores, reverse=True))

    x = np.arange(len(scores))
    x_norm = (x - x.min()) / (x.max() - x.min())
    y_norm = (scores - scores.min()) / (scores.max() - scores.min())

    # Vector from first to last point
    line_vec = np.array([x_norm[-1] - x_norm[0], y_norm[-1] - y_norm[0]])
    line_vec /= np.linalg.norm(line_vec)

    point_vecs = np.column_stack((x_norm - x_norm[0], y_norm - y_norm[0]))
    proj_lens = point_vecs @ line_vec
    proj_points = np.outer(proj_lens, line_vec)
    dists = np.linalg.norm(point_vecs - proj_points, axis=1)

    elbow_idx = np.argmax(dists)
    elbow_score = scores[elbow_idx]

    return elbow_idx, elbow_score, scores

def main():
    if len(sys.argv) < 2:
        print("Usage: detect_elbow.py bitscores.txt", file=sys.stderr)
        sys.exit(1)

    bitscore_file = sys.argv[1]
    out_plot = os.path.splitext(bitscore_file)[0] + "_elbow_hist_plot.pdf"

    # Load scores
    with open(bitscore_file) as f:
        scores = [float(line.strip()) for line in f if line.strip()]

    if len(scores) < 5:
        print("Error: too few scores to calculate elbow point.", file=sys.stderr)
        sys.exit(1)

    scores_array = np.array(scores)
    percentile_10 = np.percentile(scores_array, 10)

    elbow_idx, elbow_score, sorted_scores = detect_elbow(scores)

    # Print result
    #print(int(elbow_score))
    print(f"Elbow cutoff: {int(elbow_score)}")
    print(f"10th %tile cutoff: {int(percentile_10)}")

    # Plot sorted bit scores and histogram side-by-side
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    # Left figure: Sorted bit scores with elbow point line
    ax1.plot(sorted_scores, label="Bit scores")
    ax1.axvline(elbow_idx, color='red', linestyle='--', label=f"Elbow: {int(elbow_score)}")
    ax1.set_xlabel("Hit rank (sorted)")
    ax1.set_ylabel("Bit score")
    ax1.set_title("Bit score distribution with elbow point")
    ax1.legend()

    # Right figure: Histogram of bit scores with both cutoff lines
    ax2.hist(scores, bins=30, color='skyblue', edgecolor='black')
    ax2.axvline(elbow_score, color='red', linestyle='--', label=f"Elbow: {int(elbow_score)}")
    ax2.axvline(percentile_10, color='green', linestyle='--', label=f"10%tile: {int(percentile_10)}")
    ax2.set_xlabel("Bit score")
    ax2.set_ylabel("Frequency")
    ax2.set_title("Histogram of bit scores")
    ax2.legend()

    plt.tight_layout()
    plt.savefig(out_plot, format="pdf", bbox_inches="tight")
    print(f"Elbow + histogram plot saved to: {out_plot}", file=sys.stderr)

if __name__ == "__main__":
    main()
