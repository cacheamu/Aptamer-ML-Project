"""
free_base_count.py

Counts unpaired (free) nucleotides in aptamer secondary structures represented in
dot-bracket notation. Free bases are defined as '.' characters that occur outside
any paired region.

Input:  Text file containing one dot-bracket structure per line
Output: Text file reporting free base count for each structure

Example:
    python free_base_count.py --input secondary_structures.txt --output free_base_counts.txt
"""

import argparse


def count_free_bases(dot_bracket):
    """
    Count free bases occurring outside paired regions.
    Returns an integer.
    """
    free_count = 0
    stack = []

    for char in dot_bracket:
        if char == "(":
            stack.append("(")
        elif char == ")":
            if stack:
                stack.pop()
        elif char == ".":
            if not stack:  # count only if outside bonded region
                free_count += 1

    return free_count


def process_file(input_path, output_path):
    """
    Read input text file and write results to output text file.
    """
    with open(input_path, "r", encoding="utf-8") as f_in:
        structures = [line.strip() for line in f_in.readlines() if line.strip()]

    with open(output_path, "w", encoding="utf-8") as f_out:
        for struct in structures:
            free_bases = count_free_bases(struct)
            f_out.write(f"{struct} -> Free bases: {free_bases}\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Input text file of dot-bracket structures")
    parser.add_argument("--output", required=True, help="Output text file for free base counts")
    args = parser.parse_args()

    process_file(args.input, args.output)


if __name__ == "__main__":
    main()
