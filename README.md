# GTF_Parsing


## Description

The purpose of this script is to **parse a GTF file** using positions listed in an Excel file, in order to **annotate whether each position falls within a gene** and, if so, **which genomic feature** it overlaps.

## Input Files

- `input.gtf`  
  A GTF (Gene Transfer Format) file containing genomic features such as genes, exons, CDS, UTRs, etc.

- `input.xlsx`  
  An Excel file containing at least the following two columns:
  - `CHROM` – the chromosome name (matching the GTF)
  - `POSITION` – the genomic position of interest

## Output File

- `output.xlsx`  
  A copy of the input Excel file with **three additional columns**:
  - `GENE`: the name of the gene that overlaps the position (or `"NO"` if no gene is found)
  - `LOCATION_IN_GENE`: the specific feature within the gene, such as `exon_3`, `intron_2`, or left empty
  - `OTHER_LOCATION`: if the position is not in an exon or intron, this column will indicate special elements like:
    - `start_codon`
    - `stop_codon`
    - `five_prime_utr`
    - `three_prime_utr`
    - `CDS_2` (for coding sequence regions)

## Annotation Rules

- If the position is **outside of any gene**, `"NO"` is recorded for all annotation columns.
- If the position is **within a gene**:
  - If it overlaps an **exon**, the corresponding `exon_number` is added (e.g., `exon_1`)
  - If it overlaps a **CDS**, the corresponding `CDS_number` is added (e.g., `CDS_1`)
  - If it falls **between two exons**, it is annotated as `intron_X`, where `X` is the smaller exon number
  - If it falls in **UTRs or codons**, these are recorded in the `OTHER_LOCATION` column

## Notes

- Only the **canonical transcript** is considered (tagged as `"Ensembl_canonical"` in the GTF).
- The script assumes the GTF uses **1-based inclusive coordinates**.
- Chromosome names in both files must match (e.g., `"1"` vs `"chr1"`).

## Usage

Update the filenames at the top of the script:

```python
input_gtf = "input.gtf"
input_excel = "input.xlsx"
output_excel = "output.xlsx"
```

Then run:
```bash
python3 GTF_parser.py
```

You can install the required packages with:

```bash
pip install pandas tqdm
```
