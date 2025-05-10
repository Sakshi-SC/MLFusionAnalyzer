import os
import pandas as pd
import re
import shutil
import subprocess
import csv


# --------------------------- PART 1: Process Star-Fusion outputs --------------------------- #
# Define source and destination directories
src_dirs = ['se_output', 'pe_output']
dest_dir = 'Fusion_output'

# Create the destination directory if it doesn't exist
os.makedirs(dest_dir, exist_ok=True)

# Loop through both source directories
for src in src_dirs:
    if not os.path.exists(src):
        print(f"Source directory not found: {src}")
        continue

    for subdir in os.listdir(src):
        src_path = os.path.join(src, subdir)
        dest_path = os.path.join(dest_dir, subdir)

        if os.path.isdir(src_path):
            print(f"Moving {src_path} → {dest_path}")
            shutil.move(src_path, dest_path)
        else:
            print(f"Skipping non-directory: {src_path}")

# Change working directory to Fusion_output/
fusion_dir = os.path.join(os.getcwd(), "Fusion_output")
os.chdir(fusion_dir)

print(f"Working directory set to: {os.getcwd()}")

# Process Star-Fusion files
current_dir = os.getcwd()
dataframes = []

for dir_name in os.listdir(current_dir):
    dir_path = os.path.join(current_dir, dir_name)
    if os.path.isdir(dir_path):
        input_file_path = os.path.join(dir_path, "star-fusion.fusion_predictions.abridged.coding_effect.tsv")
        output_file_path = os.path.join(dir_path, "updated_star-fusion.fusion_predictions.abridged.coding_effect.tsv")

        if os.path.exists(input_file_path):
            df = pd.read_csv(input_file_path, sep="\t", engine="python")
            df.insert(0, "Sample_ID", dir_name)
            df.to_csv(output_file_path, sep="\t", index=False)
            print(f"Updated file created: {output_file_path}")
            dataframes.append(df)
        else:
            print(f"File not found: {input_file_path}")

if dataframes:
    merged_df = pd.concat(dataframes, ignore_index=True)
    
    # Clean and process the merged dataframe
    columns_to_remove = [
        "est_J", "est_S", "SpliceType", "LargeAnchorSupport", "LeftBreakEntropy",
        "RightBreakEntropy", "annots", "CDS_LEFT_RANGE", "CDS_RIGHT_RANGE", "FUSION_MODEL",
        "FUSION_CDS", "FUSION_TRANSL", "PFAM_LEFT", "PFAM_RIGHT"
    ]
    merged_df = merged_df.drop([col for col in columns_to_remove if col in merged_df.columns], axis=1, errors='ignore')
    
    # Process breakpoints and genes
    merged_df[['Chromosome1', 'LeftBreakpoint_Pos', 'LeftStrand']] = merged_df['LeftBreakpoint'].str.split(':', expand=True)
    merged_df[['Chromosome2', 'RightBreakpoint_Pos', 'RightStrand']] = merged_df['RightBreakpoint'].str.split(':', expand=True)
    merged_df = merged_df.drop(['LeftBreakpoint', 'RightBreakpoint'], axis=1)
    merged_df["LeftGene"] = merged_df["LeftGene"].str.extract(r"([^:^]+)$")
    merged_df["RightGene"] = merged_df["RightGene"].str.extract(r"([^:^]+)$")
    
    if '#FusionName' in merged_df.columns:
        merged_df['#FusionName'] = merged_df['#FusionName'].str.replace('--', '->', regex=False)
    
    merged_df['Splice_Site'] = merged_df['LeftBreakDinuc'] + '-' + merged_df['RightBreakDinuc']
    merged_df = merged_df.drop(['LeftBreakDinuc', 'RightBreakDinuc'], axis=1)
    
    merged_df.rename(columns={
        'LeftGene': '5_geneid',
        'RightGene': '3_geneid',
        'PROT_FUSION_TYPE': 'Splice_Pattern',
        'RightBreakpoint_Pos': 'RightBreakpoint',
        'LeftBreakpoint_Pos': 'LeftBreakpoint'
    }, inplace=True)
    
    merged_df['Total_Count_(SC+RC)'] = merged_df['JunctionReadCount'] + merged_df['SpanningFragCount']
    merged_df = merged_df.drop(['JunctionReadCount', 'SpanningFragCount'], axis=1)
    merged_df['Splice_Pattern_Class'] = merged_df['Splice_Site'].apply(lambda x: 'Canonical' if x == 'GT-AG' else 'Non-Canonical')
    
    merged_output_path = os.path.join(current_dir, "1_initial_merged_predictions.tsv")
    merged_df.to_csv(merged_output_path, sep="\t", index=False)
    print(f"Initial merged file created: {merged_output_path}")

# --------------------------- PART 2: Process GTF to extract longest transcript per gene --------------------------- #
gtf_file = "Arabidopsis_thaliana.TAIR10.51.gtf"
columns = ["seqname", "source", "feature", "start", "end", "score", "strand", "frame", "attribute"]

gtf_df = pd.read_csv(gtf_file, sep="\t", comment="#", names=columns, dtype=str)
gtf_df = gtf_df[gtf_df["feature"] == "exon"]

gtf_df["gene_id"] = gtf_df["attribute"].str.extract(r'gene_id "([^"]+)"')
gtf_df["transcript_id"] = gtf_df["attribute"].str.extract(r'transcript_id "([^"]+)"')
gtf_df["exon_number"] = gtf_df["attribute"].str.extract(r'exon_number "(\d+)"').astype(float)

max_exon_per_gene = gtf_df.groupby("gene_id")["exon_number"].max().reset_index()
result = gtf_df.merge(max_exon_per_gene, on=["gene_id", "exon_number"], how="inner")
result = result[["gene_id", "transcript_id", "exon_number"]].drop_duplicates()

output_file = "2_filtered_gtf.tsv"
result.to_csv(output_file, sep="\t", index=False, header=["Gene_ID", "Transcript_ID", "Exon_Count"])
print(f"Filtered GTF saved as {output_file}")

# --------------------------- PART 3: Fill missing transcript IDs --------------------------- #
fusion_df = pd.read_csv("1_initial_merged_predictions.tsv", sep="\t")
filtered_gtf = pd.read_csv("2_filtered_gtf.tsv", sep="\t")

gene_to_transcripts = filtered_gtf.groupby("Gene_ID")["Transcript_ID"].apply(list).to_dict()

def fill_transcript_ids(row):
    left_transcripts = gene_to_transcripts.get(row["5_geneid"], ["."]) if row["CDS_LEFT_ID"] == "." else [row["CDS_LEFT_ID"]]
    right_transcripts = gene_to_transcripts.get(row["3_geneid"], ["."]) if row["CDS_RIGHT_ID"] == "." else [row["CDS_RIGHT_ID"]]
    
    expanded_rows = []
    for left in left_transcripts:
        for right in right_transcripts:
            new_row = row.copy()
            new_row["CDS_LEFT_ID"] = left
            new_row["CDS_RIGHT_ID"] = right
            expanded_rows.append(new_row)
    return expanded_rows

expanded_data = fusion_df.apply(fill_transcript_ids, axis=1)
expanded_df = pd.DataFrame([item for sublist in expanded_data for item in sublist])

expanded_df["CDS_LEFT_ID"] = expanded_df["CDS_LEFT_ID"].str.replace("transcript:", "", regex=False)
expanded_df["CDS_RIGHT_ID"] = expanded_df["CDS_RIGHT_ID"].str.replace("transcript:", "", regex=False)

expanded_df.to_csv("3_merged_filled.tsv", sep="\t", index=False)
print("Temporary expanded DataFrame saved as '3_merged_filled.tsv'.")

# --------------------------- PART 4: Add exon information --------------------------- #
def parse_gtf(gtf_file):
    exon_map = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.strip().split("\t")
            if cols[2] == "exon":
                chrom = cols[0]
                start = int(cols[3])
                end = int(cols[4])
                strand = cols[6]
                attr = cols[8]
                transcript_match = re.search(r'transcript_id "([^"]+)"', attr)
                exon_num_match = re.search(r'exon_number "([^"]+)"', attr)
                if transcript_match and exon_num_match:
                    transcript_id = transcript_match.group(1)
                    exon_number = exon_num_match.group(1)
                    exon_map.setdefault(transcript_id, []).append({
                        'start': start,
                        'end': end,
                        'exon_number': exon_number,
                        'chrom': chrom,
                        'strand': strand
                    })
    return exon_map

print("Parsing GTF file for exon information...")
exon_data = parse_gtf(gtf_file)
print("GTF parsing completed.\n")

def find_exon(transcript_id, breakpoint, chrom):
    exons = exon_data.get(transcript_id, [])
    for exon in exons:
        if exon['chrom'] == chrom and exon['start'] <= breakpoint <= exon['end']:
            return exon['exon_number']
    return "NF"

print("Adding exon information to fusion data...")
expanded_df['Left_Exon'] = expanded_df.apply(
    lambda row: find_exon(row['CDS_LEFT_ID'], int(row['LeftBreakpoint']), str(row['Chromosome1'])),
    axis=1
)

expanded_df['Right_Exon'] = expanded_df.apply(
    lambda row: find_exon(row['CDS_RIGHT_ID'], int(row['RightBreakpoint']), str(row['Chromosome2'])),
    axis=1
)

final_output_path = os.path.join(current_dir, "4_final_merged_predictions.tsv")
expanded_df.to_csv(final_output_path, sep="\t", index=False)
print(f"Final output saved as: {final_output_path}")

if os.path.exists("3_merged_filled.tsv"):
    os.remove("3_merged_filled.tsv")
    print("Temporary file '3_merged_filled.tsv' removed.")

print("All processing completed successfully.")

# --------------------------- PART 5: Annotate Fusion Features --------------------------- #
df = pd.read_csv("4_final_merged_predictions.tsv", sep="\t")

# Initialize new columns
chromosome_feature = []
strand_classification = []
reciprocal_fusion_detected = []
fusion_dict = {}

# Step 1: Determine chromosomal/strand relationships & collect fusions
for _, row in df.iterrows():
    chromosome_feature.append("Intrachromosomal" if row['Chromosome1'] == row['Chromosome2'] else "Interchromosomal")
    strand_classification.append("Yes" if row['LeftStrand'] == row['RightStrand'] else "No")
    
    sample_id = row['Sample_ID']
    fusion_gene = row['#FusionName']
    if sample_id not in fusion_dict:
        fusion_dict[sample_id] = set()
    fusion_dict[sample_id].add(fusion_gene)

# Step 2: Check for reciprocal fusion presence
for _, row in df.iterrows():
    sample_id = row['Sample_ID']
    fusion_gene = row['#FusionName']
    reverse_fusion = "->".join(fusion_gene.split("->")[::-1])
    reciprocal_fusion_detected.append("Yes" if reverse_fusion in fusion_dict.get(sample_id, set()) else "No")

# Add results to DataFrame
df['Chromosome_Feature'] = chromosome_feature
df['Same_Strand'] = strand_classification
df['Reciprocal_Fusion'] = reciprocal_fusion_detected

df.to_csv("5_fusion_annotated.tsv", sep="\t", index=False)
print("✅ Step 1: Annotated fusion file saved as '5_fusion_annotated.tsv'")

# --------------------------- PART 6: Extract Gene Info from GTF --------------------------- #
print("⏳ Extracting gene info from GTF...")

gtf_file = "../Arabidopsis_thaliana.TAIR10.51.gtf"
output_file = "6_gene_info.tsv"

gene_data = []

with open(gtf_file, 'r') as f:
    for line in f:
        if line.startswith("#"):
            continue
        parts = line.strip().split('\t')
        if len(parts) != 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attributes = parts
        attr_dict = {}
        for attr in attributes.split(';'):
            if attr.strip():
                key_value = attr.strip().replace('"', '').split(' ')
                if len(key_value) >= 2:
                    attr_dict[key_value[0]] = ' '.join(key_value[1:])
        gene_id = attr_dict.get('gene_id', 'NA')
        gene_name = attr_dict.get('gene_name', 'NA')
        transcript_id = attr_dict.get('transcript_id', 'NA')

        if feature == "gene":
            gene_data.append([chrom, start, end, strand, gene_id, gene_name])

# Save gene-level data
pd.DataFrame(gene_data, columns=['chromosome', 'start', 'end', 'strand', 'gene_id', 'gene_name']) \
    .to_csv(output_file, sep='\t', index=False)
print(f"✅ Step 2: Gene data extracted to {output_file}")

# --------------------------- PART 7: Merge Gene Info with Fusion Data --------------------------- #
print("⏳ Merging gene info with fusion data...")

temp_df = pd.read_csv("5_fusion_annotated.tsv", sep="\t")
temp_df[["5_geneid", "3_geneid"]] = temp_df["#FusionName"].str.split("->", expand=True)

gene_info_df = pd.read_csv("6_gene_info.tsv", sep="\t")
gene_info_df.rename(columns={"gene_id": "geneid", "start": "gene_start", "end": "gene_end"}, inplace=True)
gene_info_df["gene_length"] = gene_info_df["gene_end"].astype(int) - gene_info_df["gene_start"].astype(int) + 1
gene_info_df = gene_info_df[["geneid", "gene_start", "gene_end", "gene_length"]]

# Merge to get 5' gene info
temp_df = temp_df.merge(gene_info_df, left_on="5_geneid", right_on="geneid", how="left")
temp_df.rename(columns={"gene_start": "5_gene_start", "gene_end": "5_gene_end", "gene_length": "5_gene_length"}, inplace=True)
temp_df.drop(columns=["geneid"], inplace=True)

# 3' gene
temp_df = temp_df.merge(gene_info_df, left_on="3_geneid", right_on="geneid", how="left")
temp_df.rename(columns={"gene_start": "3_gene_start", "gene_end": "3_gene_end", "gene_length": "3_gene_length"}, inplace=True)
temp_df.drop(columns=["geneid"], inplace=True)

# Convert columns to integers where applicable
int_columns = ['5_gene_start', '5_gene_end', '5_gene_length', '3_gene_start', '3_gene_end', '3_gene_length']
temp_df[int_columns] = temp_df[int_columns].astype('Int64')

temp_df.to_csv("7_fusion_with_gene_info.tsv", sep="\t", index=False)
print("✅ Step 3: Gene start/end/length merged and saved as '7_fusion_with_gene_info.tsv'")

# --------------------------- PART 8: Count Junctions Per Fusion Pair --------------------------- #
print("⏳ Calculating junction counts...")

df = pd.read_csv("7_fusion_with_gene_info.tsv", sep="\t")
df["fusion_pair"] = df["5_geneid"].astype(str) + "->" + df["3_geneid"].astype(str)

fusion_counts = df.groupby(["5_geneid", "3_geneid"]).apply(
    lambda x: x[["LeftBreakpoint", "RightBreakpoint"]].drop_duplicates().shape[0]
).reset_index(name="junction_count")

df = df.merge(fusion_counts, how="left", on=["5_geneid", "3_geneid"])
df["atlternative_junction"] = df["junction_count"].apply(lambda x: "yes" if x > 1 else "no")

# Save temporary
df.to_csv("8_fusion_with_junction_counts.tsv", sep="\t", index=False)

# Optional formatted counts
fusion_counts["formatted"] = fusion_counts["junction_count"].astype(str) + " " + fusion_counts["5_geneid"].astype(str) + "->" + fusion_counts["3_geneid"].astype(str)
fusion_counts[["formatted"]].to_csv("9_fusion_counts_formatted.tsv", index=False, header=False)

print(f"✅ Step 4: Junction counts added and saved as '8_fusion_with_junction_counts.tsv'")

# --------------------------- PART 9: Breakpoint Location Classification --------------------------- #
print("⏳ Classifying breakpoint locations relative to exons...")

def extract_exons(gtf_file):
    exon_map = {}
    with open(gtf_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            parts = line.strip().split('\t')
            if len(parts) != 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attributes = parts
            if feature != 'exon':
                continue
            attr_dict = {}
            for attr in attributes.split(';'):
                if attr.strip():
                    key_value = attr.strip().replace('"', '').split(' ')
                    if len(key_value) >= 2:
                        attr_dict[key_value[0]] = ' '.join(key_value[1:])
            transcript_id = attr_dict.get('transcript_id')
            gene_id = attr_dict.get('gene_id')
            if gene_id not in exon_map:
                exon_map[gene_id] = []
            exon_map[gene_id].append({"start": int(start), "end": int(end)})
    return exon_map

def classify_breakpoint_location(breakpoint, exons):
    for exon in exons:
        if breakpoint == exon["start"]:
            return "SE"
        elif breakpoint == exon["end"]:
            return "EE"
        elif exon["start"] < breakpoint < exon["end"]:
            return "M"
    return "O"

# Load exon map
exon_map = extract_exons(gtf_file)

# Initialize columns for output
left_exon_location = []
right_exon_location = []

for _, row in df.iterrows():
    left_bp = row['LeftBreakpoint']
    right_bp = row['RightBreakpoint']
    left_gene = row['5_geneid']
    right_gene = row['3_geneid']

    left_loc = classify_breakpoint_location(left_bp, exon_map.get(left_gene, []))
    right_loc = classify_breakpoint_location(right_bp, exon_map.get(right_gene, []))

    left_exon_location.append(left_loc)
    right_exon_location.append(right_loc)

df['5_loc'] = left_exon_location
df['3_loc'] = right_exon_location

df.to_csv("10_fusion_with_breakpoint_locations.tsv", sep='\t', index=False)

# --------------------------- PART 10: Exon Count Integration --------------------------- #
print("⏳ Incorporating exon count information...")

# Load the main data file
main_df = pd.read_csv("10_fusion_with_breakpoint_locations.tsv", sep="\t")

# Load the exon count file
exon_df = pd.read_csv("2_filtered_gtf.tsv", sep="\t", names=["Gene_ID", "Transcript_ID", "Exon_Count"])

# Aggregate exon counts for each Gene_ID
exon_counts = exon_df.groupby("Gene_ID")["Exon_Count"].unique().to_dict()

# Function to get exon count for a given gene_id
def get_exon_count(gene_id):
    return exon_counts.get(gene_id, [])

# Expand rows if exon counts are non-unique
data_expanded = []
for _, row in main_df.iterrows():
    exon5 = get_exon_count(row["5_geneid"])
    exon3 = get_exon_count(row["3_geneid"])

    for e5 in exon5:
        for e3 in exon3:
            new_row = row.copy()
            new_row["exon_count5"] = int(float(e5)) if pd.notna(e5) else ""
            new_row["exon_count3"] = int(float(e3)) if pd.notna(e3) else ""
            data_expanded.append(new_row)

# Convert back to DataFrame
expanded_df = pd.DataFrame(data_expanded)

# Save the final output
expanded_df.to_csv("11_final_fusion_annotated.tsv", sep='\t', index=False)
print("✅ FINAL STEP: Exon count added. Output saved as '11_final_fusion_annotated.tsv'")

# --------------------------- SHS Sequence Extraction and Analysis --------------------------- #
def extract_and_save(input_file):
    df = pd.read_csv(input_file, sep='\t')
    
    # Extracting relevant columns for 5' and 3' sides
    df_5 = df[['Chromosome1', 'LeftBreakpoint', '5_geneid', 'LeftStrand']].copy()
    df_5.insert(3, '.', '.')
    df_5.columns = ['Chromosome1', 'LeftBreakpoint', '5_geneid', '.', 'LeftStrand']
    df_5.to_csv('5_input', sep='\t', index=False, header=False)

    df_3 = df[['Chromosome2', 'RightBreakpoint', '3_geneid', 'RightStrand']].copy()
    df_3.insert(3, '.', '.')
    df_3.columns = ['Chromosome2', 'RightBreakpoint', '3_geneid', '.', 'RightStrand']
    df_3.to_csv('3_input', sep='\t', index=False, header=False)

def run_command(command):
    result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(f"Error: {result.stderr.decode()}")
    return result.stdout.decode()

def process_thresholds():
    thresholds = range(2, 10, 1)

    for threshold in thresholds:
        print(f"Processing threshold {threshold}...")

        run_command(f"awk -v threshold={threshold} '{{start=$2-threshold; end=$2+threshold; if (start < 0) start=0; print $1\"\\t\"start\"\\t\"end\"\\t\"$3\"\\t.\\t\"$5}}' 5_input > 5_output_threshold_{threshold}.bed")
        run_command(f"awk -v threshold={threshold} '{{start=$2-threshold; end=$2+threshold; if (start < 0) start=0; print $1\"\\t\"start\"\\t\"end\"\\t\"$3\"\\t.\\t\"$5}}' 3_input > 3_output_threshold_{threshold}.bed")

        run_command(f"bedtools getfasta -fi Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -bed 5_output_threshold_{threshold}.bed -name -s > SHS_arab_5_threshold_{threshold}.fa")
        run_command(f"bedtools getfasta -fi Arabidopsis_thaliana.TAIR10.dna.toplevel.fa -bed 3_output_threshold_{threshold}.bed -name -s > SHS_arab_3_threshold_{threshold}.fa")

        if not os.path.isfile(f"SHS_arab_5_threshold_{threshold}.fa") or not os.path.isfile(f"SHS_arab_3_threshold_{threshold}.fa"):
            print(f"Error: FASTA files for threshold {threshold} are not created.")
            continue

        with open(f"SHS_arab_3_threshold_{threshold}.fa", "r") as f_in, open(f"even_lines_3_with_genes_threshold_{threshold}.txt", "w") as f_out:
            lines = f_in.readlines()
            for i in range(0, len(lines), 2):
                gene = lines[i].strip()
                seq = lines[i + 1].strip()
                f_out.write(f"{gene} {seq}\n")

        with open(f"SHS_arab_5_threshold_{threshold}.fa", "r") as f_in, open(f"even_lines_5_with_genes_threshold_{threshold}.txt", "w") as f_out:
            lines = f_in.readlines()
            for i in range(0, len(lines), 2):
                gene = lines[i].strip()
                seq = lines[i + 1].strip()
                f_out.write(f"{gene} {seq}\n")

        if not os.path.isfile(f"even_lines_3_with_genes_threshold_{threshold}.txt") or not os.path.isfile(f"even_lines_5_with_genes_threshold_{threshold}.txt"):
            print(f"Error: Gene-sequence files for threshold {threshold} are not created.")
            continue

# Run the SHS extraction
input_file = "11_final_fusion_annotated.tsv"
extract_and_save(input_file)
process_thresholds()

# --------------------------- SHS Analysis --------------------------- #
def merge_shs_files(thresholds, output_file="12_merged_shs_all_thresholds.txt"):
    with open(output_file, "w") as out:
        out.write("5'_Gene_Name\t3'_Gene_Name\tSequence_5\tSequence_3\tThreshold\n")

        for threshold in thresholds:
            file_5 = f"even_lines_5_with_genes_threshold_{threshold}.txt"
            file_3 = f"even_lines_3_with_genes_threshold_{threshold}.txt"

            if not os.path.exists(file_5) or not os.path.exists(file_3):
                print(f"Skipping threshold {threshold}: Missing input files.")
                continue

            with open(file_5, "r") as f5, open(file_3, "r") as f3:
                lines_5 = f5.readlines()
                lines_3 = f3.readlines()

                for i in range(len(lines_5)):
                    gene_5, seq_5 = lines_5[i].strip().split(" ")
                    gene_3, seq_3 = lines_3[i].strip().split(" ")

                    out.write(f"{gene_5}\t{gene_3}\t{seq_5}\t{seq_3}\t{threshold}\n")

    print(f"✅ Merged file saved as: {output_file}")

# Define threshold range
threshold_range = range(2, 10, 1)
merge_shs_files(threshold_range)

# --------------------------- SHS Processing --------------------------- #
input_file = "12_merged_shs_all_thresholds.txt"
output_file = "13_merged_shs_with_shs_check.txt"
shs_only_output = "14_merged_shs_only_yes.txt"
highest_threshold_output = "15_merged_shs_highest_threshold.txt"

# Read the data
df = pd.read_csv(input_file, sep="\t")

# Function to extract actual coordinates
def compute_actual_coordinates(row):
    try:
        match1 = re.search(r'(\d+):(\d+)-(\d+)', row["5'_Gene_Name"])
        match2 = re.search(r'(\d+):(\d+)-(\d+)', row["3'_Gene_Name"])

        if not match1 or not match2:
            return pd.NA, pd.NA

        start1, end1 = int(match1.group(2)), int(match1.group(3))
        start2, end2 = int(match2.group(2)), int(match2.group(3))
        threshold = int(row["Threshold"])

        actual_coordinate5 = start1 + threshold
        actual_coordinate3 = end2 - threshold

        return actual_coordinate5, actual_coordinate3
    except Exception as e:
        print(f"Error processing row: {row}\nError: {e}")
        return pd.NA, pd.NA

# Apply coordinate computation
df[["Actual_Coordinate5", "Actual_Coordinate3"]] = df.apply(compute_actual_coordinates, axis=1, result_type="expand")

# Function to check SHS
def check_shs(row):
    return "YES" if row["Sequence_5"] == row["Sequence_3"] else "NO"

# Apply SHS check
df["SHS_Status"] = df.apply(check_shs, axis=1)

# Save full annotated output
df.to_csv(output_file, sep="\t", index=False)

# Filter only rows with SHS = YES
df_shs_only = df[df["SHS_Status"] == "YES"].copy()

# Extract gene names without coordinates
def extract_gene_name(full_name):
    return full_name.split("::")[0]

# Apply transformations using .loc to avoid warnings
df_shs_only.loc[:, "5'_Gene"] = df_shs_only["5'_Gene_Name"].apply(extract_gene_name)
df_shs_only.loc[:, "3'_Gene"] = df_shs_only["3'_Gene_Name"].apply(extract_gene_name)
df_shs_only.loc[:, "Threshold"] = df_shs_only["Threshold"].astype(int)

# Save SHS only file
df_shs_only.to_csv(shs_only_output, sep="\t", index=False)

# Keep rows with max threshold for each gene pair
df_highest = df_shs_only.loc[
    df_shs_only.groupby(["5'_Gene", "3'_Gene", "Actual_Coordinate5", "Actual_Coordinate3"])["Threshold"].idxmax()
]

# Save final result
df_highest.to_csv(highest_threshold_output, sep="\t", index=False)

print(f"✅ Processed file saved as: {output_file}")
print(f"✅ SHS only file saved as: {shs_only_output}")
print(f"✅ Highest threshold SHS file saved as: {highest_threshold_output}")

# --------------------------- Final SHS Annotation --------------------------- #
junction_file = "11_final_fusion_annotated.tsv"
shs_file = "15_merged_shs_highest_threshold.txt"
output_file = "tempo_fusion_data_annotated.tsv"

# Load Data
df_junction = pd.read_csv(junction_file, sep="\t")
df_shs = pd.read_csv(shs_file, sep="\t")

# Helper: Extract gene name
def extract_gene_name(full_name):
    return full_name.split("::")[0].strip(">")

# Prepare SHS Data
df_shs["5'_Gene"] = df_shs["5'_Gene_Name"].apply(extract_gene_name)
df_shs["3'_Gene"] = df_shs["3'_Gene_Name"].apply(extract_gene_name)
df_shs["Actual_Coordinate5"] = df_shs["Actual_Coordinate5"].astype(int)
df_shs["Actual_Coordinate3"] = df_shs["Actual_Coordinate3"].astype(int)
df_shs["shs_sequence"] = df_shs["Sequence_5"]
df_shs["shs_length"] = df_shs["Threshold"].astype(int) * 2

# Annotate Junctions with SHS
def find_shs(row):
    match = df_shs[
        (df_shs["5'_Gene"] == row["5_geneid"]) &
        (df_shs["3'_Gene"] == row["3_geneid"]) &
        (df_shs["Actual_Coordinate5"] == row["LeftBreakpoint"]) &
        (df_shs["Actual_Coordinate3"] == row["RightBreakpoint"])
    ]
    if not match.empty:
        return pd.Series([match.iloc[0]["shs_sequence"], match.iloc[0]["shs_length"]])
    return pd.Series(["NF", "NF"])

df_junction[["shs_sequence", "shs_length"]] = df_junction.apply(find_shs, axis=1)

# Modify Specific SHS Sequences
def modify_shs(row):
    if str(row["shs_length"]) == "4" and "AGGT" in str(row["shs_sequence"]):
        row["shs_sequence"] = row["shs_sequence"].replace("AGGT", "NF")
        row["shs_length"] = "NF"
    return row

df_junction = df_junction.apply(modify_shs, axis=1)

# Save Output
df_junction.to_csv(output_file, sep="\t", index=False)
print(f"✅ File processed successfully! Saved as {output_file}")

#------total_mapped_reads-------


# Load the fusion data
fusion_df = pd.read_csv("tempo_fusion_data_annotated.tsv", sep="\t")

# Load the reads count data
reads_df = pd.read_csv("reads.tsv", sep="\t")

# Merge based on sample ID
merged_df = fusion_df.merge(reads_df, left_on="Sample_ID", right_on="sample", how="left")

# Drop the extra 'sample' column (optional)
merged_df.drop(columns=["sample"], inplace=True)

# Save the updated file
merged_df.to_csv("tempo_fusion_data_annotated_with_reads.tsv", sep="\t", index=False)

print("Updated file saved as: tempo_fusion_data_annotated_with_reads.tsv")


—-

input_file = 'tempo_fusion_data_annotated.tsv'
output_file = 'fusion_data_with_mapped_reads.tsv'

with open(input_file, 'r') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter='\t')
    fieldnames = reader.fieldnames + ['Total_Mapped_Reads']
    writer = csv.DictWriter(outfile, fieldnames=fieldnames, delimiter='\t')
    writer.writeheader()

    for row in reader:
        ffpm = row['FFPM']
        total_count = row['Total_Count_(SC+RC)']
        
        if ffpm not in ['NA', 'NF', ''] and total_count not in ['NA', 'NF', '']:
            try:
                ffpm_value = float(ffpm)
                count_value = int(total_count)
                if ffpm_value > 0:
                    total_reads = (count_value * 1_000_000) / ffpm_value
                    row['Total_Mapped_Reads'] = int(total_reads)
                else:
                    row['Total_Mapped_Reads'] = 'NF'
            except ValueError:
                row['Total_Mapped_Reads'] = 'NF'
        else:
            row['Total_Mapped_Reads'] = 'NF'
        
        writer.writerow(row)

print(f"✅ Done. Output saved to: {output_file}")



# --------------------------- GO Analysis --------------------------- #

# Step 1: Load fusion data
df = pd.read_csv("fusion_data_with_mapped_reads.tsv", sep="\t")

# Step 2: Extract gene IDs and save to separate files
five_prime_ids = df["5_geneid"].dropna().unique().tolist()
three_prime_ids = df["3_geneid"].dropna().unique().tolist()

# Save to files

with open("gene_ids_5prime.txt", "w") as f:
    f.write("\n".join(five_prime_ids))
print("✅ 5' gene IDs saved to 'gene_ids_5prime.txt'")


with open("gene_ids_3prime.txt", "w") as f:
    f.write("\n".join(three_prime_ids))
print("✅ 3' gene IDs saved to 'gene_ids_3prime.txt'")




 
import requests
import json

# Step 1: Read gene IDs from the input files
def read_gene_file(filename):
    with open(filename, 'r') as f:
        gene_ids = [line.strip() for line in f if line.strip()]
    return gene_ids

# Load your gene lists
gene_ids_3prime = read_gene_file('gene_ids_3prime.txt')
gene_ids_5prime = read_gene_file('gene_ids_5prime.txt')

# Step 2: Function to query g:Profiler
def run_gprofiler_query(query_data, organism='osativa'):
    url = 'https://biit.cs.ut.ee/gprofiler/api/gost/profile/'
    payload = {
        'organism': organism,
        'query': query_data,
        'sources': ['GO:BP'],  #'GO:MF','GO:CC
        'user_threshold': 0.05,
        'significance_threshold_method': 'g_SCS',
        'all_results': False,
        'ordered': False,
        'no_iea': True,
        'no_evidences': False,  # Keep evidence codes for now
        'output': 'json'
    }
    response = requests.post(url, json=payload, headers={'User-Agent': 'PythonRequest'})
    if response.status_code != 200:
        print(f"Error: API request failed with status code {response.status_code}")
        print(response.text)
        return None
    return response.json()['result']

# Step 3: Process and save results
def save_enrichment_results(results, query_name, input_gene_ids, output_filename):
    with open(output_filename, 'w') as f:
        # Write header
        f.write("Query\tTerm_ID\tTerm_Name\tP_Value\tIntersection_Genes\tIntersection_Size\tTerm_Size\tQuery_Size\tSource\n")
        
        # Process each result
        for result in results:
            if result['query'] == query_name and result['significant']:
                # Get intersecting genes by mapping non-empty intersections
                intersection_data = result.get('intersections', [])
                intersecting_genes = []

                if intersection_data and len(intersection_data) == len(input_gene_ids):
                    for i, intersection in enumerate(intersection_data):
                        if intersection:  # Non-empty list (e.g., ["IDA"]) indicates this gene intersects
                            intersecting_genes.append(input_gene_ids[i])
                else:
                    print(f"WARNING: Intersections length ({len(intersection_data)}) doesn’t match input size ({len(input_gene_ids)})")
                
                # Debug: Print what we found
                print(f"DEBUG: Term {result['native']} ({result['name']}):")
                print(f"  Intersection_Size: {result['intersection_size']}")
                print(f"  Intersecting Genes: {intersecting_genes}")
                
                # Join gene IDs into a string
                intersection_gene_str = ",".join(intersecting_genes) if intersecting_genes else "N/A"
                
                # Write the row
                f.write(f"{result['query']}\t"
                        f"{result['native']}\t"
                        f"{result['name']}\t"
                        f"{result['p_value']:.2e}\t"
                        f"{intersection_gene_str}\t"
                        f"{result['intersection_size']}\t"
                        f"{result['term_size']}\t"
                        f"{result['query_size']}\t"
                        f"{result['source']}\n")

# Step 4: Run queries separately
print("Running query for 5prime_genes...")
results_5prime = run_gprofiler_query(gene_ids_5prime)
if results_5prime:
    save_enrichment_results(results_5prime, 'query_1', gene_ids_5prime, 'go_enrichment_5prime.txt')

print("Running query for 3prime_genes...")
results_3prime = run_gprofiler_query(gene_ids_3prime)
if results_3prime:
    save_enrichment_results(results_3prime, 'query_1', gene_ids_3prime, 'go_enrichment_3prime.txt')

print("Enrichment analysis complete. Results saved to 'go_enrichment_3prime.txt' and 'go_enrichment_5prime.txt'.")


# Step 1: Read go_enrichment_*.txt files and extract GO term mappings
bp_files = ["go_enrichment_3prime.txt", "go_enrichment_5prime.txt"]
go_dict = {}

for file in bp_files:
    # Read tab-separated file
    df = pd.read_csv(file, sep="\t")
    for _, row in df.iterrows():
        go_term_name = row["Term_Name"]  # Use Term_Name from g:Profiler output
        genes = row["Intersection_Genes"].split(",")  # Use Intersection_Genes from g:Profiler output
        for gene in genes:
            if gene in go_dict:
                go_dict[gene].append(go_term_name)
            else:
                go_dict[gene] = [go_term_name]

# Step 2: Read the fusion gene file (fusion_data_annotated.tsv)
fusion_df = pd.read_csv("temp_fusion_data_annotated.tsv", sep="\t")  # Updated to tab-separated TSV

# Step 3: Map GO terms to 5_geneid and 3_geneid
def get_go_terms(gene_id):
    return "|".join(go_dict.get(gene_id, ["No_GO_Term"]))

fusion_df["5_gene_go_term"] = fusion_df["5_geneid"].apply(get_go_terms)
fusion_df["3_gene_go_term"] = fusion_df["3_geneid"].apply(get_go_terms)

# Step 4: Save the updated file
fusion_df.to_csv("fusion_annotated.tsv", index=False, sep="\t")  # Updated output name

print("Processing complete. File saved as “fusion_data_annotated.tsv")

