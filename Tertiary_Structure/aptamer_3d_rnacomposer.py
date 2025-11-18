"""
aptamer_3d_rnacomposer.py

Sends RNA secondary structure (dot-bracket) and sequence to RNAComposer API
to generate a 3D PDB model of an RNA aptamer.

Requirements:
    pip install requests

Note:
    Output is a .pdb file containing tertiary structure coordinates.
"""

import requests


def submit_to_rnacomposer(sequence, structure, output_file="aptamer_3d.pdb"):
    url = "http://rnacomposer.cs.put.poznan.pl/rest"
    payload = {
        "sequence": sequence,
        "structure": structure
    }
    response = requests.post(url + "/generate", json=payload)

    if response.status_code != 200:
        raise Exception("Error submitting job to RNAComposer")

    job_id = response.text.strip()
    result_url = f"{url}/result/{job_id}"

    # Fetch the output
    model = requests.get(result_url)
    with open(output_file, "w") as f:
        f.write(model.text)

    print(f"3D model saved to {output_file}")


if __name__ == "__main__":
    # Replace with your aptamer
    sequence = "AGCUAGCUAGGCGUUA"
    structure = "((((...(((....)))))))"   # dot-bracket from RNAfold

    submit_to_rnacomposer(sequence, structure, output_file="aptamer_3d.pdb")
