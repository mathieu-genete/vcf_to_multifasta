# VCF to Multifasta Script (Version 2.0.0)

This script converts VCF files to multifasta format, incorporating population information, ploidy levels, and reference sequences for greater flexibility and accuracy. It can handle compressed VCF files.

---

## Features

- Support for VCF files (both `.vcf` and `.vcf.gz` formats).
- Population-aware multifasta generation.
- Reference sequence integration and validation.
- Filtering options using contig/chromosome lists.
- Random allele selection for heterozygous genotypes.
- Customizable ploidy settings.

---

## Installation

1. Clone or download the repository containing this script.
2. Ensure Python 3 and the required dependencies are installed:
   - `Biopython`
   - `PyVCF`

To install dependencies, run:

```bash
pip install biopython
pip install pyvcf
```

---

## Usage

Run the script with the following command-line options:

| **Argument**       | **Description**                                   | **Required**      | **Default**      |
|----------------------|--------------------------------------------------|-------------------|------------------|
| `-v, --vcffile`      | Input VCF file (.vcf or .vcf.gz).                | Required          | N/A              |
| `-o, --outfolder`    | Output folder for fasta files.                   | Required          | N/A              |
| `-r, --reffasta`     | Fasta reference file.                            | Optional          | None             |
| `-c, --contigslist`  | Text file with contig or chromosome names.       | Optional          | All contigs.     |
| `-n, --popname`      | File with individual and population names.       | Optional          | None             |
| `-s, --randomSeed`   | Random seed for reproducibility.                 | Optional          | Current time.    |
| `-p, --ploidy`       | Ploidy of samples.                               | Optional          | `2`.             |
| `-nw, --hidewarnings`| hide warnings messages                           | Optional          | False            |

---

## Example Command Line

Convert a VCF file to multifasta format using ploidy 2 and reference sequences:

```bash
python vcf_to_multifasta_v2.py -v data/sample.vcf.gz -o output/ \
                               -r data/reference.fasta \
                               -c data/contigs.txt \
                               -n data/populations.txt \
                               -p 2 -s 12345
```

If a reference fasta is not provided, the script will generate sequences with "N" for missing positions.

---

## Output

The script generates multifasta files in the specified output folder. Each file corresponds to a contig/chromosome and contains sequences for individuals, formatted as:

```
>{gene}|{pop}|{indiv}|allele{anbr}
ACTGACTGACTG...
```
