"""
calculate_zscores.py

Computes Z-scores for MFE (Minimum Free Energy) values obtained from
secondary structure prediction.

Input:  CSV file containing an MFE column (last column recommended)
Output: CSV file with Aptamer, MFE, and Z-score

Example:
    python calculate_zscores.py --input Final_Results.csv --output Final_Zscores.csv
"""

import argparse
import csv
import statistics


def read_mfe_values(filepath):
    """
    Read a CSV file and extract aptamer names and MFE values.
    Returns two lists: aptamer_names, mfe_values.
    """
    aptamers = []
    mfes = []

    with open(filepath, "r", encoding="utf-8") as f:
        next(f)  # skip header
        for line in f:
            parts = line.strip().split(",")
            if len(parts) >= 3:
                aptamers.append(parts[0])
                try:
                    mfes.append(float(parts[-1]))
                except ValueError:
                    mfes.append(0.0)

    return aptamers, mfes


def compute_zscores(values):
    """
    Compute Z-scores for a list of numeric values.
    """
    mean = statistics.mean(values)
    stdev = statistics.pstdev(values)  # population standard deviation

    if stdev == 0:
        return [0.0 for _ in values]

    return [(v - mean) / stdev for v in values]


def write_output(filepath, aptamers, mfes, zscores):
    """
    Write results to CSV.
    """
    with open(filepath, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Aptamer", "MFE", "Z-score"])
        for a, mfe, z in zip(aptamers, mfes, zscores):
            writer.writerow([a, f"{mfe:.2f}", f"{z:.4f}"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input CSV file with MFE column")
    parser.add_argument("--output", required=True, help="Output CSV file to save Z-scores")
    args = parser.parse_args()

    aptamers, mfes = read_mfe_values(args.input)
    zscores = compute_zscores(mfes)
    write_output(args.output, aptamers, mfes, zscores)


if __name__ == "__main__":
    main()
