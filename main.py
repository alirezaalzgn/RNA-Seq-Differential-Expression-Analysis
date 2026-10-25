# Imports
# ------------------------------------------
import pandas as pd
from pathlib import Path
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler


# ------------------------------------------
# Data Reading
# ------------------------------------------

BASE_DIR = Path(__file__).resolve().parent
DATA_DIR = BASE_DIR / "Data" / "Raw"

counts_file = DATA_DIR / "GSE316391_RNA_counts_Met_plus_GTEX.csv.gz"

counts = pd.read_csv(
    counts_file,
    low_memory=False
)


# ------------------------------------------
# Metadata Reading
# ------------------------------------------

metadata_file = (
    DATA_DIR /
    "GSE316391_BCMetastasis_All_Samples_and_Clinical_Data.xlsx"
)

excel_file = pd.ExcelFile(metadata_file)


# ------------------------------------------
# RNA-seq Metadata
# ------------------------------------------

rna_metadata = pd.read_excel(
    metadata_file,
    sheet_name="RNA-seq",
    header=None
)

# Find the actual header row
header_row = 0

for i in range(min(5, len(rna_metadata))):

    if rna_metadata.iloc[i].notna().sum() >= 3:
        header_row = i
        break


# Set the correct header
rna_metadata.columns = rna_metadata.iloc[header_row]

# Remove the header row from the data
rna_metadata = rna_metadata.iloc[header_row + 1:]

# Remove completely empty rows
rna_metadata = rna_metadata.dropna(how="all")

# Reset index
rna_metadata = rna_metadata.reset_index(drop=True)


# ------------------------------------------
# Sample Mapping Preparation
# ------------------------------------------

# Get RNA-seq sample names from metadata
metadata_samples = (
    rna_metadata["RNA-seq name"]
    .astype(str)
    .tolist()
)


# Get sample columns from count matrix
count_sample_columns = [
    col
    for col in counts.columns
    if col not in [
        "gene_name",
        "Geneid",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length",
        "Description.x",
        "Description.y",
        "Description"
    ]
]


print("\n===== Sample Mapping Check =====")

print(
    "Number of metadata samples:",
    len(metadata_samples)
)

print(
    "Number of count matrix samples:",
    len(count_sample_columns)
)


# ------------------------------------------
# Build Sample Mapping
# ------------------------------------------

def convert_rna_name(rna_name):
    """
    Convert RNA-seq metadata name
    to count matrix sample name.
    """

    name = str(rna_name)

    # Handle known sample naming discrepancy
    if name == "R49LIV_2":
        return "R49_4-Liver_Tumor"

    # Replace separators used in metadata
    name = name.replace("-", "_")

    # Extract patient/sample part and tissue
    if "BR_" in name:

        patient, replicate = name.split("BR_")
        tissue = "Breast"

    elif "LIV_" in name:

        patient, replicate = name.split("LIV_")
        tissue = "Liver"

    elif "LUN_" in name:

        patient, replicate = name.split("LUN_")
        tissue = "Lung"

    else:

        return None

    return f"{patient}_{replicate}-{tissue}_Tumor"


# Create mapped sample names
rna_metadata["Count_Sample"] = (
    rna_metadata["RNA-seq name"]
    .apply(convert_rna_name)
)


# ------------------------------------------
# Validate Sample Mapping
# ------------------------------------------

count_sample_set = set(count_sample_columns)

rna_metadata["Mapping_OK"] = (
    rna_metadata["Count_Sample"]
    .isin(count_sample_set)
)


# ------------------------------------------
# Mapping Summary
# ------------------------------------------

print("\n===== Mapping Summary =====")

print(
    "Total metadata samples:",
    len(rna_metadata)
)

print(
    "Successfully mapped:",
    rna_metadata["Mapping_OK"].sum()
)

print(
    "Unmapped:",
    (~rna_metadata["Mapping_OK"]).sum()
)


# ------------------------------------------
# Show Unmapped Samples
# ------------------------------------------

unmapped_samples = rna_metadata[
    ~rna_metadata["Mapping_OK"]
]

if len(unmapped_samples) > 0:

    print("\n===== Unmapped Samples =====")

    print(
        unmapped_samples[
            [
                "RNA-seq name",
                "Patient #",
                "Tissue Site",
                "Count_Sample"
            ]
        ].to_string(index=False)
    )


# ------------------------------------------
# Build Analysis Metadata
# ------------------------------------------

analysis_metadata = rna_metadata[
    [
        "RNA-seq name",
        "Patient #",
        "Tissue Site",
        "External ID",
        "Container ID",
        "Internal Case ID",
        "Count_Sample",
        "Mapping_OK",
        "Histology",
        "Sequencing Instrument",
        "Sequencing Details"
    ]
].copy()


# Simplify tissue labels
analysis_metadata["Tissue"] = (
    analysis_metadata["Tissue Site"]
    .replace(
        {
            "Liver/Bile Duct": "Liver",
            "Breast": "Breast",
            "Lung": "Lung"
        }
    )
)


# ------------------------------------------
# Analysis Metadata Check
# ------------------------------------------

print("\n===== Analysis Metadata =====")

print(
    "Shape:",
    analysis_metadata.shape
)

print("\nTissue counts:")

print(
    analysis_metadata["Tissue"].value_counts()
)


# ------------------------------------------
# Clinical Metadata
# ------------------------------------------

clinical_metadata = pd.read_excel(
    metadata_file,
    sheet_name="Clinical Data"
).dropna(how="all")


# ------------------------------------------
# Prepare Merge Keys
# ------------------------------------------

analysis_metadata["Patient #"] = pd.to_numeric(
    analysis_metadata["Patient #"],
    errors="coerce"
)

clinical_metadata["Patient #"] = pd.to_numeric(
    clinical_metadata["Patient #"],
    errors="coerce"
)


# ------------------------------------------
# Prepare Clinical Metadata for Merge
# ------------------------------------------

clinical_metadata = clinical_metadata.rename(
    columns={
        "Study ID": "External ID"
    }
)


# ------------------------------------------
# Define Merge Keys
# ------------------------------------------

merge_keys = [
    "Patient #",
    "Tissue Site",
    "External ID",
    "Container ID",
    "Internal Case ID"
]


# ------------------------------------------
# Check Clinical Metadata Uniqueness
# ------------------------------------------

print("\n===== Clinical Metadata Key Check =====")

duplicate_count = (
    clinical_metadata
    .duplicated(subset=merge_keys)
    .sum()
)

print(
    "Duplicate clinical rows:",
    duplicate_count
)


# ------------------------------------------
# Merge RNA-seq and Clinical Metadata
# ------------------------------------------

merged_metadata = analysis_metadata.merge(
    clinical_metadata,
    on=merge_keys,
    how="left",
    suffixes=("_RNA", "_Clinical"),
    indicator=True
)


# ------------------------------------------
# Merge Validation
# ------------------------------------------

print("\n===== Clinical Merge Validation =====")

print(
    "RNA-seq samples:",
    len(analysis_metadata)
)

print(
    "Merged rows:",
    len(merged_metadata)
)

print("\nMerge status:")

print(
    merged_metadata["_merge"].value_counts()
)


# ------------------------------------------
# Show Unmatched Samples
# ------------------------------------------

unmatched = merged_metadata[
    merged_metadata["_merge"] != "both"
]

if len(unmatched) == 0:

    print(
        "\nAll RNA-seq samples matched successfully."
    )

else:

    print("\n===== Unmatched Samples =====")

    print(
        unmatched[
            [
                "RNA-seq name",
                "Patient #",
                "Tissue Site",
                "External ID",
                "Container ID",
                "Internal Case ID"
            ]
        ].to_string(index=False)
    )


# ------------------------------------------
# Build Study Count Matrix
# ------------------------------------------

# Keep only successfully mapped study samples
study_samples = analysis_metadata[
    analysis_metadata["Mapping_OK"] == True
]["Count_Sample"].tolist()


# Select gene annotation columns
gene_columns = [
    col
    for col in [
        "gene_name",
        "Geneid",
        "Chr",
        "Start",
        "End",
        "Strand",
        "Length",
        "Description.x",
        "Description.y",
        "Description"
    ]
    if col in counts.columns
]


# Keep gene information + study samples
study_count_matrix = counts[
    gene_columns + study_samples
].copy()


print("\n===== Study Count Matrix =====")

print(
    "Shape:",
    study_count_matrix.shape
)

print(
    "Number of study samples:",
    len(study_samples)
)


# ------------------------------------------
# Library Size QC
# ------------------------------------------

# Calculate total counts for each study sample
library_sizes = (
    study_count_matrix[study_samples]
    .sum(axis=0)
    .sort_values(ascending=False)
)


print("\n===== Library Size Summary =====")

print(
    library_sizes.describe()
)


# ------------------------------------------
# Gene Filtering QC
# ------------------------------------------

# Extract raw count matrix
count_data = study_count_matrix[study_samples].copy()


# Count genes expressed above different thresholds
genes_at_least_1 = (count_data >= 1).sum(axis=1)
genes_at_least_10 = (count_data >= 10).sum(axis=1)
genes_at_least_20 = (count_data >= 20).sum(axis=1)


print("\n===== Gene Filtering QC =====")

print(
    "Genes with count >= 1 in at least 1 sample:",
    (genes_at_least_1 >= 1).sum()
)

print(
    "Genes with count >= 10 in at least 1 sample:",
    (genes_at_least_10 >= 1).sum()
)

print(
    "Genes with count >= 20 in at least 1 sample:",
    (genes_at_least_20 >= 1).sum()
)


# Genes expressed in at least 50% of samples
min_samples = len(study_samples) // 2

genes_50_percent = (
    (count_data >= 10).sum(axis=1) >= min_samples
)


print(
    f"Genes with count >= 10 in at least {min_samples} samples:",
    genes_50_percent.sum()
)


# ------------------------------------------
# PCA - Sample Level QC
# ------------------------------------------

# Raw count matrix
count_data = study_count_matrix[study_samples].copy()


# Keep genes expressed in at least 50% of samples
min_samples = len(study_samples) // 2

filter_mask = (
    (count_data >= 10).sum(axis=1) >= min_samples
)

filtered_counts = count_data.loc[filter_mask].copy()


print("\n===== PCA QC =====")

print(
    "Genes before filtering:",
    count_data.shape[0]
)

print(
    "Genes after filtering:",
    filtered_counts.shape[0]
)


# Log-transform counts for exploratory PCA
log_counts = np.log2(filtered_counts + 1)


# Transpose: samples × genes
pca_input = log_counts.T


# Standardize genes
scaler = StandardScaler()

pca_scaled = scaler.fit_transform(pca_input)


# PCA
pca = PCA(n_components=2)

pca_result = pca.fit_transform(pca_scaled)


# Explained variance
explained_variance = pca.explained_variance_ratio_


print(
    "Explained variance (PC1, PC2):",
    explained_variance[0],
    explained_variance[1]
)


# Create PCA dataframe
pca_df = pd.DataFrame(
    pca_result,
    columns=["PC1", "PC2"]
)

pca_df["Sample"] = study_samples


# Add metadata
sample_metadata = (
    analysis_metadata[
        ["Count_Sample", "Patient #", "Tissue"]
    ]
    .drop_duplicates()
)

pca_df = pca_df.merge(
    sample_metadata,
    left_on="Sample",
    right_on="Count_Sample",
    how="left"
)


# Plot PCA
plt.figure(figsize=(10, 7))

for tissue in pca_df["Tissue"].unique():

    subset = pca_df[
        pca_df["Tissue"] == tissue
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=tissue,
        s=80
    )


plt.xlabel(
    f"PC1 ({explained_variance[0] * 100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({explained_variance[1] * 100:.1f}%)"
)

plt.title("PCA of RNA-seq Samples")

plt.legend()

plt.tight_layout()

plt.show()


# ------------------------------------------
# PCA - Normalized Expression QC
# ------------------------------------------

# Raw count matrix
count_data = study_count_matrix[study_samples].copy()


# Keep genes expressed in at least 50% of samples
min_samples = len(study_samples) // 2

filter_mask = (
    (count_data >= 10).sum(axis=1) >= min_samples
)

filtered_counts = count_data.loc[filter_mask].copy()


# Calculate library sizes
library_sizes = filtered_counts.sum(axis=0)


# Convert counts to CPM
cpm = filtered_counts.div(
    library_sizes,
    axis=1
) * 1_000_000


# Log-transform normalized expression
log_cpm = np.log2(cpm + 1)


# Transpose: samples × genes
pca_input = log_cpm.T


# PCA without per-gene standardization
pca = PCA(n_components=2)

pca_result = pca.fit_transform(pca_input)


# Explained variance
explained_variance = pca.explained_variance_ratio_


print("\n===== Normalized PCA QC =====")

print(
    "Genes used:",
    filtered_counts.shape[0]
)

print(
    "Explained variance (PC1, PC2):",
    explained_variance[0],
    explained_variance[1]
)


# Create PCA dataframe
pca_df = pd.DataFrame(
    pca_result,
    columns=["PC1", "PC2"]
)

pca_df["Sample"] = study_samples


# Add metadata
sample_metadata = (
    analysis_metadata[
        ["Count_Sample", "Patient #", "Tissue"]
    ]
    .drop_duplicates()
)


pca_df = pca_df.merge(
    sample_metadata,
    left_on="Sample",
    right_on="Count_Sample",
    how="left"
)


# Plot PCA
plt.figure(figsize=(10, 7))

for tissue in pca_df["Tissue"].unique():

    subset = pca_df[
        pca_df["Tissue"] == tissue
    ]

    plt.scatter(
        subset["PC1"],
        subset["PC2"],
        label=tissue,
        s=80
    )


plt.xlabel(
    f"PC1 ({explained_variance[0] * 100:.1f}%)"
)

plt.ylabel(
    f"PC2 ({explained_variance[1] * 100:.1f}%)"
)

plt.title(
    "PCA of RNA-seq Samples - Log2 CPM"
)

plt.legend()

plt.tight_layout()

plt.show()


# ------------------------------------------
# Sample-to-Sample Correlation
# ------------------------------------------

# Use the normalized log-CPM matrix from PCA
correlation_matrix = log_cpm.corr()


print("\n===== Sample-to-Sample Correlation =====")

print(
    "Correlation matrix shape:",
    correlation_matrix.shape
)

print(
    "Min / Max correlation:",
    correlation_matrix.values[
        ~np.eye(correlation_matrix.shape[0], dtype=bool)
    ].min(),
    correlation_matrix.values[
        ~np.eye(correlation_matrix.shape[0], dtype=bool)
    ].max()
)


# ------------------------------------------
# Replicate Structure Check
# ------------------------------------------

# Create replicate groups
replicate_groups = (
    analysis_metadata
    .groupby(["Patient #", "Tissue"])
    ["Count_Sample"]
    .apply(list)
    .reset_index(name="Samples")
)


# Keep groups with more than one sample
replicate_groups = replicate_groups[
    replicate_groups["Samples"].apply(len) > 1
].copy()


print("\n===== Replicate Groups =====")

for _, row in replicate_groups.iterrows():

    patient = row["Patient #"]
    tissue = row["Tissue"]
    samples = row["Samples"]

    print(
        f"Patient {patient} - {tissue}:",
        samples
    )


# ------------------------------------------
# Build DESeq2 Input Data
# ------------------------------------------

print("\n===== Building DESeq2 Input =====")


# ------------------------------------------
# Select paired Liver/Lung patients
# ------------------------------------------

# Patients with both Liver and Lung samples
paired_patients = [1, 3, 4, 5, 6, 7]


# Keep only Liver and Lung samples
de_metadata = analysis_metadata[
    analysis_metadata["Patient #"].isin(paired_patients)
    & analysis_metadata["Tissue"].isin(["Liver", "Lung"])
].copy()


print(
    "Patients selected for paired analysis:",
    sorted(de_metadata["Patient #"].unique())
)

print(
    "Tissue counts:",
    de_metadata["Tissue"].value_counts().to_dict()
)


# ------------------------------------------
# Create biological sample groups
# ------------------------------------------

de_metadata["Biological_Sample"] = (
    "P"
    + de_metadata["Patient #"].astype(int).astype(str)
    + "_"
    + de_metadata["Tissue"]
)


# ------------------------------------------
# Collapse technical replicates
# ------------------------------------------

collapsed_counts = {}

for sample_name, group in de_metadata.groupby(
    "Biological_Sample"
):

    replicate_columns = group["Count_Sample"].tolist()

    collapsed_counts[sample_name] = (
        count_data[replicate_columns]
        .sum(axis=1)
    )


# Convert to DataFrame
de_counts = pd.DataFrame(
    collapsed_counts,
    index=count_data.index
)


# ------------------------------------------
# Add gene identifiers
# ------------------------------------------

de_counts.insert(
    0,
    "Geneid",
    study_count_matrix["Geneid"].values
)


de_counts.insert(
    1,
    "gene_name",
    study_count_matrix["gene_name"].values
)


# ------------------------------------------
# Build DE metadata
# ------------------------------------------

de_sample_metadata = (
    de_metadata
    .groupby("Biological_Sample")
    .first()
    .reset_index()
)


de_sample_metadata = de_sample_metadata[
    [
        "Biological_Sample",
        "Patient #",
        "Tissue"
    ]
].copy()


de_sample_metadata = de_sample_metadata.rename(
    columns={
        "Biological_Sample": "Sample",
        "Patient #": "Patient",
        "Tissue": "Tissue"
    }
)


# Set categorical variables
de_sample_metadata["Patient"] = (
    de_sample_metadata["Patient"]
    .astype(int)
    .astype(str)
)

de_sample_metadata["Tissue"] = (
    de_sample_metadata["Tissue"]
    .astype("category")
)


# ------------------------------------------
# Reorder samples
# ------------------------------------------

sample_order = []

for patient in paired_patients:

    for tissue in ["Liver", "Lung"]:

        sample_name = (
            f"P{patient}_{tissue}"
        )

        sample_order.append(sample_name)


de_counts = de_counts[
    ["Geneid", "gene_name"] + sample_order
]


de_sample_metadata = (
    de_sample_metadata
    .set_index("Sample")
    .loc[sample_order]
    .reset_index()
)


# ------------------------------------------
# Save DESeq2 input files
# ------------------------------------------

processed_dir = Path(
    "Data/Processed"
)

processed_dir.mkdir(
    parents=True,
    exist_ok=True
)


de_counts_path = (
    processed_dir /
    "DE_counts.csv"
)


de_metadata_path = (
    processed_dir /
    "DE_metadata.csv"
)


de_counts.to_csv(
    de_counts_path,
    index=False
)


de_sample_metadata.to_csv(
    de_metadata_path,
    index=False
)


print("\n===== DESeq2 Input Saved =====")

print(
    "DE count matrix shape:",
    de_counts.shape
)

print(
    "DE metadata shape:",
    de_sample_metadata.shape
)



# Within-Group Replicate Correlation
# ------------------------------------------

print("\n===== Within-Group Replicate Correlation =====")

for _, row in replicate_groups.iterrows():
    patient = row["Patient #"]
    tissue = row["Tissue"]
    samples = row["Samples"]

    group_counts = count_data[samples]

    # Use a separate variable to preserve the full sample correlation matrix
    group_correlation_matrix = group_counts.corr(method="pearson")

    correlations = []

    for i in range(len(samples)):
        for j in range(i + 1, len(samples)):
            correlations.append(
                group_correlation_matrix.iloc[i, j]
            )

    print(
        f"Patient {patient} - {tissue}:",
        [round(x, 4) for x in correlations]
    )

# ------------------------------------------
# Correlation Heatmap
# ------------------------------------------

plt.figure(figsize=(14, 12))

sns.heatmap(
    correlation_matrix,
    cmap="viridis",
    square=True,
    xticklabels=True,
    yticklabels=True
)

plt.title("Sample-to-Sample Correlation")

plt.xticks(
    rotation=90,
    fontsize=8
)

plt.yticks(
    rotation=0,
    fontsize=8
)

plt.tight_layout()

plt.show()


# ------------------------------------------
# Hierarchical Clustering
# ------------------------------------------

sns.clustermap(
    correlation_matrix,
    cmap="viridis",
    figsize=(14, 14),
    xticklabels=True,
    yticklabels=True
)

plt.show()


# ------------------------------------------
# Final DESeq2 Input QC
# ------------------------------------------

print("\n===== Final DESeq2 Input QC =====")

count_columns = sample_order

print(
    "Missing values (counts / metadata):",
    de_counts.isna().sum().sum(),
    de_sample_metadata.isna().sum().sum()
)

print(
    "Negative count values:",
    (de_counts[count_columns] < 0).sum().sum()
)

print(
    "Non-integer count values:",
    de_counts[count_columns]
    .map(lambda x: x % 1 != 0)
    .sum().sum()
)

print(
    "Duplicate Gene IDs:",
    de_counts["Geneid"].duplicated().sum()
)

print(
    "Missing Gene IDs:",
    de_counts["Geneid"].isna().sum()
)


# ------------------------------------------
# Resolve Duplicate Gene IDs
# ------------------------------------------

before_rows = len(de_counts)

de_counts = de_counts.drop_duplicates(
    subset="Geneid",
    keep="first"
).copy()

after_rows = len(de_counts)

print("\n===== Resolving Duplicate Gene IDs =====")

print(
    "Rows before / after removing duplicates:",
    before_rows,
    after_rows
)

print(
    "Duplicate rows removed:",
    before_rows - after_rows
)


# ------------------------------------------
# Save Clean DESeq2 Input
# ------------------------------------------

de_counts.to_csv(
    de_counts_path,
    index=False
)

de_sample_metadata.to_csv(
    de_metadata_path,
    index=False
)

print("\n===== Clean DESeq2 Input Saved =====")

print("DE counts:", de_counts_path)
print("DE metadata:", de_metadata_path)

