PCR Candidate Scanner + Primer3 Pipeline

Files:

main.py — runs the complete workflow
500bp_scanner.py — scans and ranks overlapping 500 bp regions
primer3_design.py — designs and ranks Primer3 primer pairs
requirements.txt — Python dependency list

Installation:
py -m pip install -r requirements.txt

Run:
python main.py
Paste a FASTA or DNA sequence, then press Enter on a blank line.

Output: 
The program can save:
scanner_candidates.csv
ranked_primer_pairs.csv
ranked_candidate_regions.fasta
Notes

The filename 500bp_scanner.py begins with a number, so main.py loads itusing Python's importlib rather than a normal import statement.