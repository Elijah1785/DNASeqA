# DNASeqA

**Automated DNA Candidate Region Screening and PCR Primer Design**

DNASeqA is a Python-based bioinformatics pipeline that identifies promising regions within DNA sequences and designs PCR primer pairs using Primer3.

The project combines sequence-quality assessment, automated primer design, and custom ranking methods into a single workflow to help streamline preliminary PCR assay development.

## Features

- **DNA region screening:** Scans DNA sequences using overlapping 500 bp windows.
- **Sequence-quality evaluation:** Assesses GC content, sequence complexity, homopolymer runs, repetitive sequences, and ambiguous bases.
- **Primer3 integration:** Generates candidate PCR primer pairs for selected regions.
- **Custom ranking:** Evaluates candidate regions and primer pairs using a combination of sequence-quality and primer-design metrics.
- **Interactive results:** Allows users to inspect ranked results and individual primer pairs.
- **Exportable results:** Saves candidate-region and primer-design results in CSV and FASTA formats.

## Project Structure

```text
DNASeqA/
├── src/
│   ├── main.py
│   ├── region_scanner.py
│   └── primer3_design.py
├── tests/
│   ├── conftest.py
│   ├── test_region_scanner.py
│   └── test_primer3_design.py
├── README.md
├── README_TESTS.md
├── requirements.txt
└── requirements-dev.txt
```

## Installation

### Requirements

- Python
- pip
- Dependencies listed in `requirements.txt`

Clone the repository:

```bash
git clone https://github.com/Elijah1785/DNASeqA.git
cd DNASeqA
```

Install dependencies:

```bash
python -m pip install -r requirements.txt
```

## Usage

Run the main program from the project root:

```bash
python src/main.py
```

When prompted, paste a DNA sequence or FASTA-formatted sequence. Press Enter on a blank line to finish entering the sequence.

DNASeqA will scan the input sequence, evaluate candidate regions, and attempt to design primer pairs using Primer3.

### Interactive Commands

After processing, users can explore results with the following commands:

| Command | Description |
|---|---|
| `list` | Display ranked results |
| `regions` | Display scanner candidates |
| `show 1` | Display detailed information for result 1 |
| `save` | Export results |
| `quit` | Exit the program |

## Output Files

DNASeqA can export:

| File | Contents |
|---|---|
| `scanner_candidates.csv` | Candidate regions and associated screening metrics |
| `ranked_primer_pairs.csv` | Ranked primer pairs and design metrics |
| `ranked_candidate_regions.fasta` | Candidate DNA region sequences in FASTA format |

These files can be used for additional analysis, documentation, or downstream research workflows.

## How It Works

DNASeqA follows a multi-stage computational workflow:

1. **Sequence input:** Accepts DNA sequence data.
2. **Region scanning:** Examines overlapping 500 bp windows.
3. **Quality assessment:** Evaluates GC content, entropy, repeats, homopolymers, and ambiguous bases.
4. **Candidate ranking:** Prioritizes regions based on sequence-quality metrics.
5. **Primer design:** Uses Primer3 to generate candidate PCR primer pairs.
6. **Combined scoring:** Ranks results using region quality and primer quality.
7. **Result inspection and export:** Provides interactive output and downloadable data files.

The combined ranking uses a weighted scoring approach to prioritize candidate regions with favorable primer-design characteristics.

## Testing

The project includes automated tests for sequence screening and primer-design functionality.

Install development dependencies:

```bash
python -m pip install -r requirements-dev.txt
```

Run the tests:

```bash
python -m pytest -v
```

Additional testing information is available in `README_TESTS.md`.

## Limitations

DNASeqA is intended for preliminary computational screening and primer design.

- Primer rankings are heuristic and do not guarantee successful PCR amplification.
- The pipeline does not currently perform genome-wide primer specificity verification.
- Candidate primer pairs may require additional evaluation for off-target binding.
- Experimental validation is necessary before primers are used in laboratory applications.

DNASeqA uses Primer3 as its underlying primer-design engine. Its primary contribution is the integration of region screening, custom ranking, and primer design into a unified workflow.

## Future Development

Potential future improvements include:

- BLAST-based primer specificity assessment.
- Expanded command-line options.
- Support for batch sequence processing.
- Additional reporting and visualization capabilities.
- Improved installation and distribution options.

## Author

Developed by Elijah Kleiman.

## Project Status

**Version 0.1.0 — Initial public development release**

DNASeqA is under active development, and additional features and improvements may be introduced in future versions.
