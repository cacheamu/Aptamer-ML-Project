"""
fetch_properties_from_pubchem.py

Retrieves physicochemical properties from PubChem using CID
(Compound ID retrieved from: fetch_cids_pubchempy.py)

Input:  Text file containing one CID per line
Output: CSV file containing selected compound properties

Properties retrieved:
    - Molecular weight
    - XLogP
    - Hydrogen bond donor count
    - Hydrogen bond acceptor count
    - Rotatable bond count

Example:
    python fetch_properties_from_pubchem.py --input cids.txt --output properties.csv
"""

import argparse
import csv
import pubchempy as pcp


def fetch_properties(cid):
    """
    Retrieve key physicochemical properties from PubChem.
    Returns a dictionary of property values.
    """
    try:
        compound = pcp.Compound.from_cid(cid)
        return {
            "CID": cid,
            "Molecular weight": compound.molecular_weight,
            "XLogP": compound.xlogp,
            "HBD": compound.h_bond_donor_count,
            "HBA": compound.h_bond_acceptor_count,
            "Rotatable bonds": compound.rotatable_bond_count
        }
    except Exception:
        return None


def process_file(input_path, output_path):
    """
    Read CIDs and write compound properties to CSV.
    """
    with open(input_path, "r", encoding="utf-8") as f:
        cids = [line.strip() for line in f.readlines() if line.strip()]

    results = []
    for cid in cids:
        try:
            numeric_cid = int(cid)
        except ValueError:
            continue

        props = fetch_properties(numeric_cid)
        if props:
            results.append(props)

    with open(output_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f,
            fieldnames=[
                "CID",
                "Molecular weight",
                "XLogP",
                "HBD",
                "HBA",
                "Rotatable bonds"
            ]
        )
        writer.writeheader()
        writer.writerows(results)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input text file containing CIDs")
    parser.add_argument("--output", required=True, help="Output CSV file for properties")
    args = parser.parse_args()

    process_file(args.input, args.output)


if __name__ == "__main__":
    main()
