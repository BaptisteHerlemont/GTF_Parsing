import pandas as pd
import re
from collections import defaultdict
from tqdm import tqdm

input_gtf = "input.gtf"
input_excel = "input.xlsx"
output_excel = "output.xlsx"

df = pd.read_excel(input_excel)

def parse_attributes(attr_str):
    attrs = {}
    for part in attr_str.strip().split(";"):
        if part.strip() == '':
            continue
        key, value = part.strip().split(' ', 1)
        attrs[key] = value.replace('"', '').strip()
    return attrs

gene_infos = defaultdict(lambda: defaultdict(dict))

with open(input_gtf) as gtf:
    for line in gtf:
        if line.startswith("#"):
            continue
        cols = line.strip().split('\t')
        if len(cols) != 9:
            continue
        chrom, source, feature, start, end, _, strand, _, attrs = cols
        attrs = parse_attributes(attrs)
        if feature == "transcript":
            if "tag" in attrs and "Ensembl_canonical" in attrs["tag"]:
                t_id = attrs["transcript_id"]
                g_id = attrs["gene_id"]
                g_name = attrs.get("gene_name", g_id)
                gene_infos[chrom][g_id]["transcript_id"] = t_id
                gene_infos[chrom][g_id]["gene_name"] = g_name
                gene_infos[chrom][g_id]["gene_start"] = int(start)
                gene_infos[chrom][g_id]["gene_end"] = int(end)
                gene_infos[chrom][g_id]["strand"] = strand
                gene_infos[chrom][g_id]["exons"] = []
                gene_infos[chrom][g_id]["cds"] = []
                gene_infos[chrom][g_id]["five_prime_utr"] = []
                gene_infos[chrom][g_id]["three_prime_utr"] = []
                gene_infos[chrom][g_id]["start_codon"] = []
                gene_infos[chrom][g_id]["stop_codon"] = []
        elif feature in ["exon", "CDS", "five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"]:
            t_id = attrs.get("transcript_id", None)
            g_id = attrs.get("gene_id", None)
            if not t_id or not g_id:
                continue
            if g_id in gene_infos[chrom] and gene_infos[chrom][g_id].get("transcript_id", None) == t_id:
                if feature == "exon":
                    gene_infos[chrom][g_id]["exons"].append((int(start), int(end), attrs.get("exon_number")))
                elif feature == "CDS":
                    gene_infos[chrom][g_id]["cds"].append((int(start), int(end), attrs.get("exon_number")))
                else:
                    gene_infos[chrom][g_id][feature].append((int(start), int(end)))

def annotate_variant(chrom, pos):
    if chrom not in gene_infos:
        return "NO", "", ""
    for g_id, info in gene_infos[chrom].items():
        if info["gene_start"] <= pos <= info["gene_end"]:
            gene_name = info["gene_name"]
            strand = info["strand"]
            for feature in ["five_prime_utr", "three_prime_utr", "start_codon", "stop_codon"]:
                for (s, e) in info[feature]:
                    if s <= pos <= e:
                        return gene_name, "", feature
            for (s, e, exon_num) in info["exons"]:
                if s <= pos <= e:
                    in_cds = False
                    for (cs, ce, cds_num) in info["cds"]:
                        if cs <= pos <= ce:
                            in_cds = True
                            break
                    if in_cds:
                        return gene_name, f"exon_{exon_num}", f"CDS_{exon_num}"
                    else:
                        return gene_name, f"exon_{exon_num}", ""
            for (cs, ce, cds_num) in info["cds"]:
                if cs <= pos <= ce:
                    return gene_name, "", f"CDS_{cds_num}"
            exons = info["exons"]
            exons_sorted = sorted(exons, key=lambda x: int(x[0]))
            for i in range(len(exons_sorted)-1):
                left = exons_sorted[i]
                right = exons_sorted[i+1]
                left_end = left[1]
                right_start = right[0]
                if (left_end < pos < right_start) or (right_start < pos < left_end):
                    exon_nums = [int(left[2]), int(right[2])]
                    intron_num = min(exon_nums)
                    return gene_name, f"intron_{intron_num}", ""
            return gene_name, "", ""
    return "NO", "", ""

GENE_list = []
LOCATION_list = []
OTHER_LOCATION_list = []

for idx, row in tqdm(df.iterrows(), total=len(df)):
    chrom = str(row["CHROM"])
    pos = int(row["POSITION"])
    gene, location, other = annotate_variant(chrom, pos)
    GENE_list.append(gene)
    LOCATION_list.append(location)
    OTHER_LOCATION_list.append(other)

df["GENE"] = GENE_list
df["LOCATION_IN_GENE"] = LOCATION_list
df["OTHER_LOCATION"] = OTHER_LOCATION_list

df.to_excel(output_excel, index=False)

