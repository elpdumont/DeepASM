# File management
import json
import logging
import os
import random
import sys

import dask.dataframe as dd
import google.cloud.logging

# Python packages for data, stats
import numpy as np
import pandas as pd
import yaml
from google.cloud import bigquery, storage
from google.cloud.logging.handlers import CloudLoggingHandler
from sklearn.neighbors import KernelDensity

from gcp_utils import create_df_from_json_for_index_file

# Initialize random seed (for selecting the reads used in the matrix)
random_seed = 42
random.seed(random_seed)

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()

# Set up logging (using GCP)
log_name = "procee_json"

# Create a handler for Google Cloud Logging.
gcloud_logging_client = google.cloud.logging.Client()
gcloud_logging_handler = CloudLoggingHandler(gcloud_logging_client, name=log_name)

# Create a stream handler to log messages to the console.
stream_handler = logging.StreamHandler()
stream_handler.setLevel(logging.WARNING)

# Now create a logger and add the handlers:
logger = logging.getLogger(log_name)
logger.setLevel(logging.DEBUG)
logger.addHandler(gcloud_logging_handler)
logger.addHandler(stream_handler)

# Import all other variables from the config file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
bucket_name = config["GCP"]["BUCKET_NAME"]

# Variables to handle
vars_to_remove = config["VARS_TO_REMOVE"]
dic_vars_to_keep = config["VARS_TO_KEEP"]
vars_to_keep = list(dic_vars_to_keep.keys())

# vars_to_normalize = config["VARS_TO_NORMALIZE"]
categorical_vars = config["CATEGORICAL_VARS"]
categorical_vars_ohe = config["CAT_VARS_OHE"]

# Genomic variables
genomic_length = config["GENOMICS"]["GENOMIC_LENGTH"]
min_nb_reads_in_sequence = config["GENOMICS"]["MIN_NB_READS_IN_SEQUENCE"]
min_fraction_of_nb_cpg_in_read = config["GENOMICS"]["MIN_FRACTION_OF_NB_CPG_IN_READ"]

# Feature prep (Kernel functions)
kernel_fm_nb_values = config["FEATURE_PREP"]["KERNEL_FM_NB_VALUES"]
kernel_fm_bandwidth = config["FEATURE_PREP"]["KERNEL_FM_BANDWIDTH"]
kernel_cov_nb_max = config["FEATURE_PREP"]["KERNEL_COV_NB_MAX"]
kernel_cov_nb_step = config["FEATURE_PREP"]["KERNEL_COV_NB_STEP"]
kernel_cov_bandwidth = config["FEATURE_PREP"]["KERNEL_COV_BANDWIDTH"]
kernel_type = config["FEATURE_PREP"]["KERNEL_TYPE"]


# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
ml_dataset_id = os.getenv("ML_DATASET_ID")
raw_data_bucket_folder = os.getenv("CLOUDASM_DATASET_ID")
nb_files_per_task = int(os.getenv("NB_FILES_PER_TASK", 0))


def create_schema_fields(variables_dict):
    """Create schema fields from a dictionary of variable names and types."""
    schema_fields = [
        bigquery.SchemaField(name, type_, mode="REQUIRED")
        for name, type_ in variables_dict.items()
    ]
    logger.info(f"Schema fields created: {schema_fields}")
    return schema_fields


def add_record_field(schema_fields, record_name, fields, mode="REPEATED"):
    """Add a record field to an existing schema."""
    record_field = bigquery.SchemaField(record_name, "RECORD", mode=mode, fields=fields)
    return schema_fields + [record_field]


def upload_dataframe(bq_client, dataframe, table_id, schema):
    """Upload a dataframe to BigQuery with the specified schema."""
    logger.info(f"Uploading dataframe to {table_id}")
    job_config = bigquery.LoadJobConfig(schema=schema)
    job = bq_client.load_table_from_dataframe(
        dataframe, table_id, job_config=job_config
    )
    result = job.result()  # Wait for the job to complete
    logger.info(f"Load job result for {table_id}: {result}")


def compute_counts_and_percentages(column):
    """
    Computes the counts and percentages of zeros and ones in a column.

    Args:
    - column (pd.Series): A pandas Series containing the data.

    Returns:
    - counts (dict): A dictionary with counts of zeros and ones.
    - percentages (dict): A dictionary with percentages of zeros and ones.
    """

    column = column.astype(int)
    # Count the absolute number of zeros and ones
    counts = column.value_counts().reindex([0, 1], fill_value=0)

    # Calculate the percentages
    percentages = (counts / len(column) * 100).round(2)

    return counts.to_dict(), percentages.to_dict()


def filter_chromosomes(df):
    """
    Ffilters out rows where the 'chr' column is 'X' or 'Y'.

    Parameters:
    - df (pd.DataFrame): The input DataFrame from which the 'sample' column will be removed and rows with 'chr' values of 'X' or 'Y' will be filtered out.

    Returns:
    - pd.DataFrame: The modified DataFrame without the 'sample' column and without rows where 'chr' is 'X' or 'Y'.
    """

    df = df[~df["chr"].isin(["X", "Y"])].copy()
    df["chr"] = df["chr"].astype(int)
    return df


def calculate_cpg_distances(cpg_positions):
    """
    Calculates the consecutive distances between CpG positions in a list.

    Parameters:
    - cpg_positions (list of int): A list of integer positions representing CpG sites.

    Returns:
    - list of int: A list of distances between consecutive CpG positions.
    """
    # Initialize an empty list to store distances between consecutive CpG positions
    distances = []

    # Iterate over pairs of consecutive CpG positions
    for current_pos, next_pos in zip(cpg_positions, cpg_positions[1:]):
        # Calculate the distance between the current position and the next, and append to the list
        distances.append(next_pos - current_pos)

    return distances


def generate_kernel_values(start, stop, step=1, reshape=True):
    """
    Generates an array of values to be used for kernel density estimation.

    Parameters:
    - start: The starting value of the range.
    - stop: The end value of the range.
    - step: The increment between each value in the range.
    - reshape: Boolean indicating whether to reshape the array for sklearn's KernelDensity.

    Returns:
    - A numpy array of values, optionally reshaped for KernelDensity input.
    """
    values = np.arange(start, stop + step, step)
    if reshape:
        values = values.reshape(-1, 1)
    return values


def estimate_kernel_density(x, values_for_kernel, bandwidth, kernel=kernel_type):
    """
    Estimates the kernel density for a given sample and set of values.

    Parameters:
    - x: The input sample data for kernel density estimation.
    - values_for_kernel: The values over which the kernel density is estimated.
    - bandwidth: The bandwidth of the kernel.
    - kernel: The type of kernel to use (default: 'gaussian').

    Returns:
    - A numpy array of the estimated densities, rounded to 4 decimal places.
    """
    sample = np.reshape(x, (-1, 1))
    kernel_model = KernelDensity(bandwidth=bandwidth, kernel=kernel)
    kernel_model.fit(sample)
    log_densities = kernel_model.score_samples(values_for_kernel)
    probabilities = np.exp(log_densities)
    return np.round(probabilities, 4)


def convert_arrays_to_kernel_densities(df, column_name, values_for_kernel, bandwidth):
    """
    Calculates the mean, standard deviation, and kernel density estimate for a given column in a DataFrame.

    Parameters:
    - df: DataFrame containing the target column.
    - column_name: Name of the column to process.
    - values_for_kernel: Pre-computed values for kernel density estimation.
    - bandwidth: The bandwidth to use for kernel density estimation.
    """
    # Calculate mean and standard deviation
    std_name = f"std_{column_name}"
    mean_name = f"mean_{column_name}"
    df[std_name] = df[column_name].apply(
        lambda x: np.round(np.std(x, ddof=1), 4)
    )  # ddof=1 for sample standard deviation
    df[mean_name] = df[column_name].apply(lambda x: np.round(np.mean(x), 4))

    # Kernel density estimates
    kernel_name = column_name + "_kd"
    df[kernel_name] = df[column_name].apply(
        lambda x: estimate_kernel_density(x, values_for_kernel, bandwidth)
    )


def expand_array_elements(df, column_name):
    """
    Expands the elements of an array from a specified column in a DataFrame into separate columns.
    Each new column will have a suffix `_k` where `k` is the index of the element in the array.

    Parameters:
    - df (pd.DataFrame): The DataFrame to be modified.
    - column_name (str): The name of the column in `df` where each element is an array.
    """

    # Check if the column exists in the DataFrame
    if column_name not in df.columns:
        print(f"Column {column_name} not found in DataFrame.")
        return

    # Determine the maximum length of the arrays in the specified column
    max_length = df[column_name].apply(len).max()

    # Generate new columns for each element in the array
    for i in range(max_length):
        # Create a new column name with the suffix `_k`
        new_column_name = f"{column_name}_{i}"

        # Extract the i-th element from each array in the specified column
        df[new_column_name] = df[column_name].apply(
            lambda x, i=i: x[i] if i < len(x) else None
        )

    # Note: This function modifies the DataFrame in place and does not return anything.


def generate_nucleotide_sequence_of_cpg_fm(row, genomic_length):
    """
    Generates a sequence representation of a genomic region, detailing the presence of CpG sites and their
    fractional methylation status. It creates a list of dictionaries, each representing a genomic position.
    The presence of a CpG site at a position is marked with 1, otherwise 0, and the fractional methylation value
    is specified for positions with CpG sites. For positions without CpG sites, the fractional methylation is
    represented as None, indicating the absence of data.

    Parameters:
    - row (pd.Series): A Pandas Series object containing 'cpg_pos' and 'cpg_fm' keys/columns. 'cpg_pos' should be
                       a list of 1-indexed positions indicating where CpG sites are located, and 'cpg_fm' should
                       be a corresponding list of fractional methylation values for these positions.
    - genomic_length (int): The total length of the genomic region to be analyzed, indicating the number of positions.

    Returns:
    - list of dict: A list of dictionaries where each dictionary represents a genomic position. Each dictionary
                    has the following keys:
                    'pos' (int): The genomic position, 1-indexed.
                    'cpg' (int): Indicates the presence of a CpG site at the position (1 for presence, 0 for absence).
                    'cpg_fm' (float or None): The fractional methylation value at positions with a CpG site; None
                                              for positions without CpG sites.
    """

    # Initialize the list to hold dictionaries for each genomic position
    result = [{"pos": k + 1, "cpg": 0, "cpg_fm": None} for k in range(genomic_length)]

    # Update dictionaries in the result list for positions with CpGs
    for pos, fm in zip(row["cpg_pos"], row["cpg_fm"]):
        # Adjust for 0-based indexing and update only the necessary keys
        index = pos - 1  # Adjusting for 0-based indexing
        result[index]["cpg"] = 1
        result[index]["cpg_fm"] = np.round(fm, 3)

    return result


def generate_nucleotide_sequence_of_cpg_and_methyl_for_read(
    read_dictionary, genomic_length, min_nb_cpg_in_read
):
    """
    Generates a methylation profile for a given read.

    The profile is a one-dimensional numpy array representing each position within
    the genomic length. Each position in the array can have a value of:
    - 0 if there is no CpG at that position,
    - 1 if there is a CpG with no methylation, and
    - 2 if there is a CpG with methylation.

    Parameters:
    - read_dictionary (dict): A dictionary containing 'pos_array' with 1-indexed positions of CpGs,
                              'meth_array' indicating methylation status (0 or 1) for each CpG,
                              'read_fm' indicating forward or reverse read, and 'nb_cpg' the number of CpGs.
    - genomic_length (int): The length of the genome or region being considered.
    - min_nb_cpg_in_read (int): Minimum number of CpGs required in a read for it to be considered.

    Returns:
    - dict: A dictionary containing 'read_fm' and 'cpg_methylation_profile' if nb_cpg meets or exceeds the minimum.
            Returns None if the condition is not met.
    """

    pos_array = (
        np.array(read_dictionary["pos_array"], dtype=int) - 1
    )  # Convert to 0-indexed
    meth_array = np.array(read_dictionary["meth_array"], dtype=int)

    # Initialize the result array with zeros
    result = np.zeros(genomic_length, dtype=int)

    # Set the values based on CpG presence and methylation status
    result[pos_array] = meth_array + 1

    return (
        {
            "read_fm": read_dictionary["read_fm"],
            "cpg_methylation_profile": result,
        }
        if int(read_dictionary["nb_cpg"]) >= min_nb_cpg_in_read
        else None
    )


def generate_sequence_cpg_cov_and_methyl_over_reads(
    row_df, genomic_length, nb_reads_in_sequence, min_fraction_of_nb_cpg
):
    """
    Generates a genomic-length wide list of dictionaries.
    Each dictionary is a dictionary with 2 keys: nucleotide position, and reads.
    Each "reads" is also dictionary with 2 keys: read number and CpG status:
       - 0 if there is no CpG at that position,
       - 1 if there is a CpG with no methylation, and
       - 2 if there is a CpG with methylation.

    Reads are ordered by fractional methylation (FM) in descending order..
    If the minimum coverage is not met, an empty list is returned.

    Parameters:
    - row_df (pd.Series): A row from a DataFrame, representing data for a specific genomic window.
    - genomic_interval (int): The length of the genomic window to be considered.
    - nb_reads_in_sequence (int): The minimum number of reads required to generate the matrix.
    - min_fraction_of_nb_cpg (float): The minimum fraction of CpGs required to select a read to be considered.

    Returns:
    - list: A list of dictionaries with 2 keys: nucleotide position, and reads.
            Each "reads" is a list of dictionaries with 2 keys: read number and CpG status:
               - 0 if there is no CpG at that position,
               - 1 if there is a CpG with no methylation, and
               - 2 if there is a CpG with methylation.
            Reads are ordered by fractional methylation (FM) in descending order across the genomic sequence.
            If the minimum nb of reads in sequence is not met, an empty list is returned.

    """

    genomic_array = row_df["genomic_picture"]
    nb_cpg_found_in_region = int(row_df["nb_cpg_found"])
    min_cpg_to_keep_in_read = int(nb_cpg_found_in_region * min_fraction_of_nb_cpg)

    # Generate the 1D CpG methylation profile for each read
    profiles = [
        generate_nucleotide_sequence_of_cpg_and_methyl_for_read(
            read_dict,
            genomic_length,
            min_cpg_to_keep_in_read,
        )
        for read_dict in genomic_array
    ]

    profiles = [profile for profile in profiles if profile]

    if len(profiles) >= nb_reads_in_sequence:

        # Randomly sample profiles if more are available than the minimum required
        sampled_profiles = (
            random.sample(profiles, nb_reads_in_sequence)
            if len(profiles) > nb_reads_in_sequence
            else profiles
        )

        # Sort the sampled profiles by read FM, descending
        sorted_profiles = sorted(
            sampled_profiles, key=lambda d: d["read_fm"], reverse=True
        )

        # Extract the CpG methylation profiles
        methylation_matrix = np.array(
            [profile["cpg_methylation_profile"] for profile in sorted_profiles]
        )

        # Transpose the matrix to match the required dimensions to be genomic_length x nb_reads_in_sequence
        methylation_matrix = np.transpose(methylation_matrix)

        pos_reads_array = []

        for pos in range(genomic_length):
            reads_info = []
            for read_nb in range(nb_reads_in_sequence):
                cpg_state = methylation_matrix[pos, read_nb]
                reads_info.append({"read_nb": read_nb + 1, "cpg_state": cpg_state})

            pos_reads_array.append({"pos": pos + 1, "reads": reads_info})

    else:
        pos_reads_array = []

    return pos_reads_array


# Define main script
def main():

    logger.info(f"Config file : {config}")

    # Store the JSON file into a dataframe
    df_raw, file_name = create_df_from_json_for_index_file(
        bucket_name, raw_data_bucket_folder, BATCH_TASK_INDEX, nb_files_per_task
    )
    logger.info(f"File names: {file_name}")
    logger.info(f"Number of rows in raw dataframe: {len(df_raw)}")

    logger.info("Create kernel functions")
    values_for_kernel_cov = generate_kernel_values(
        0, kernel_cov_nb_max, kernel_cov_nb_step
    )

    # Kernel function for fractional methylation ('read_fm', 'cpg_fm')
    values_for_kernel_fm = generate_kernel_values(0, 1, step=1 / kernel_fm_nb_values)

    dic_kernel = {
        "read_fm": {"values": values_for_kernel_fm, "bandwidth": kernel_fm_bandwidth},
        "cpg_fm": {"values": values_for_kernel_fm, "bandwidth": kernel_fm_bandwidth},
        "cpg_cov": {"values": values_for_kernel_cov, "bandwidth": kernel_cov_bandwidth},
        "cpg_dist": {
            "values": values_for_kernel_cov,
            "bandwidth": kernel_cov_bandwidth,
        },
    }

    logger.info("Remove unwanted chromosomes")
    df_filtered = filter_chromosomes(df_raw)
    logger.info(f"Number of rows of df_filtered: {len(df_filtered)}")
    logger.info(compute_counts_and_percentages(df_filtered["asm_snp"]))

    logger.info("Calculate consecutive distances between CpGs")
    df_filtered["cpg_dist"] = df_filtered["cpg_pos"].apply(calculate_cpg_distances)

    logger.info("Convert specific arrays into kernel densities")
    for col in dic_kernel.keys():
        if col in df_filtered.columns:
            convert_arrays_to_kernel_densities(
                df_filtered,
                col,
                dic_kernel[col]["values"],
                dic_kernel[col]["bandwidth"],
            )

    logger.info(
        "Create a column for each kernel density estimate (right now they are stored in arrays)"
    )
    for col in dic_kernel.keys():
        expand_array_elements(df_filtered, col + "_kd")

    logger.info(
        "Create a nucleotide sequence of CpG fractional methylation over the genomic length. Each nucleotide comes with 2 infos: presence of CpG, presence of methylation"
    )
    df_filtered["sequence_cpg_fm"] = df_filtered.apply(
        lambda x: generate_nucleotide_sequence_of_cpg_fm(x, genomic_length), axis=1
    )

    logger.info(
        "Create a sequence of nucleotide over the genomic length. Each nucleotide is represented by a vector of the depth (Number of reads required). Each element represents the presence of a CpG and its methylation status (0: no CpG, 1: Cpg non methylated, 2: methylated CpG)"
    )

    ddf = dd.from_pandas(df_filtered, npartitions=10)
    result = ddf.apply(
        lambda x: generate_sequence_cpg_cov_and_methyl_over_reads(
            x, genomic_length, min_nb_reads_in_sequence, min_fraction_of_nb_cpg_in_read
        ),
        axis=1,
        meta=("x", "object"),
    )

    df_filtered["sequence_cpg_cov_and_methyl"] = result.compute()

    initial_row_count = len(df_filtered)
    # Remove rows where the specified column contains an empty list
    df_filtered = df_filtered[
        df_filtered["sequence_cpg_cov_and_methyl"].apply(lambda x: len(x) > 0)
    ]

    # Count the number of rows after removing rows with empty lists
    final_row_count = len(df_filtered)

    # Calculate the number of rows removed
    rows_removed = initial_row_count - final_row_count

    logger.info(
        f"There are {rows_removed} rows without sequence_cpg_cov_and_methyl data from the dataset out of {initial_row_count} rows"
    )
    percentage_removed = (rows_removed / initial_row_count) * 100
    logger.info(f"{percentage_removed:.2f}% of the rows were removed")

    logger.info(compute_counts_and_percentages(df_filtered["asm_snp"]))

    # Store the different datasets into a hash table.
    dic_data = {}

    df_filtered = df_filtered.drop(
        columns=vars_to_remove, axis=1, errors="ignore"
    ).copy(deep=True)

    logger.info("One-hot encode categoricals variables that are not binary")
    dummies_list = []
    for var in categorical_vars_ohe:
        logger.info(f"One hot for variable {var}")
        dummies = pd.get_dummies(df_filtered[var], prefix=var, dtype=int)
        dummies_list.append(dummies)
    dic_data["clean"] = pd.concat(
        [df_filtered, pd.concat(dummies_list, axis=1)], axis=1
    )

    logger.info("Enforcing data types for integer variables")
    for var in [
        "chr",
        "region_inf",
        "region_sup",
        "region_nb_cpg",
        "nb_reads",
        "nb_cpg_found",
        "dnase",
        "encode_ChiP_V2",
        "tf_motifs",
    ]:
        dic_data["clean"][var] = dic_data["clean"][var].astype(pd.Int64Dtype())

    logger.info("Creating the dataset with a methylation sequence")
    dic_data["sequence_cpg_fm"] = dic_data["clean"][
        vars_to_keep + ["sequence_cpg_fm"]
    ].copy(deep=True)

    logger.info("Creating the  dataset with the methylation matrix")
    dic_data["sequence_cpg_cov_and_methyl"] = dic_data["clean"][
        vars_to_keep + ["sequence_cpg_cov_and_methyl"]
    ].copy(deep=True)

    logger.info("Creating the tabular dataset")
    dic_data["tabular"] = (
        dic_data["clean"]
        .drop(columns=["sequence_cpg_fm", "sequence_cpg_cov_and_methyl"], axis=1)
        .copy(deep=True)
    )

    logger.info(f"All variables in tabular dataset: {dic_data['tabular'].columns}")

    logger.info("Uploading the 3 datasets to BigQuery")
    # Define schema fields based on a dictionary of variables
    schema_fields = create_schema_fields(dic_vars_to_keep)

    # Define additional record fields
    record_fields_sequence_cpg_fm = [
        bigquery.SchemaField("pos", "INTEGER", mode="REQUIRED"),
        bigquery.SchemaField("cpg", "INTEGER", mode="REQUIRED"),
        bigquery.SchemaField("cpg_fm", "FLOAT", mode="NULLABLE"),
    ]

    record_fields_sequence_cpg_cov_and_methyl = [
        bigquery.SchemaField("pos", "INTEGER", mode="REQUIRED"),
        bigquery.SchemaField(
            "reads",
            "RECORD",
            mode="REPEATED",
            fields=[
                bigquery.SchemaField("read_nb", "INTEGER", mode="REQUIRED"),
                bigquery.SchemaField("cpg_state", "INTEGER", mode="REQUIRED"),
            ],
        ),
    ]

    # Add record fields to the base schema
    schema_sequence_cpg_fm = add_record_field(
        schema_fields, "sequence_cpg_fm", record_fields_sequence_cpg_fm
    )

    schema_sequence_cpg_cov_and_methyl = add_record_field(
        schema_fields,
        "sequence_cpg_cov_and_methyl",
        record_fields_sequence_cpg_cov_and_methyl,
    )

    upload_dataframe(
        bq_client,
        dic_data["sequence_cpg_fm"],
        f"{ml_dataset_id}.sequence_cpg_fm",
        schema_sequence_cpg_fm,
    )

    upload_dataframe(
        bq_client,
        dic_data["sequence_cpg_cov_and_methyl"],
        f"{ml_dataset_id}.sequence_cpg_cov_and_methyl",
        schema_sequence_cpg_cov_and_methyl,
    )

    # For tabular data with autodetection
    logger.info("Uploading tabular data with schema autodetection")
    job_config = bigquery.LoadJobConfig(autodetect=True)
    job = bq_client.load_table_from_dataframe(
        dic_data["tabular"], f"{ml_dataset_id}.tabular", job_config=job_config
    )
    logger.info(job.result())

    logger.info("SCRIPT COMPLETE")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
