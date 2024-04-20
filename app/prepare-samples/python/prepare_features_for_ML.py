# File management
import json
import logging
import os
import random
import sys

# Python packages for data, stats
import numpy as np
import pandas as pd
import yaml
from gcp import create_df_from_json_for_index_file, upload_dataframe_to_bq
from google.cloud import bigquery, storage
from sklearn.neighbors import KernelDensity

# import dask.dataframe as dd


# Initialize random seed (for selecting the reads used in the matrix)
random_seed = 42
random.seed(random_seed)

# Create a handler for Google Cloud Logging.
logging.basicConfig(level=logging.INFO)

# Import all other variables from the config file
with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
project = config["GCP"]["PROJECT"]
bucket = config["GCP"]["BUCKET"]

# Variables to handle
vars_to_remove = config["VARS_TO_REMOVE"]
dic_vars_to_keep = config["VARS_TO_KEEP"]
vars_to_keep = list(dic_vars_to_keep.keys())

categorical_vars = config["CATEGORICAL_VARS"]
categorical_vars_ohe = config["CAT_VARS_OHE"]

# Genomic variables
genomic_length = config["GENOMICS"]["GENOMIC_LENGTH"]
min_nb_reads_in_sequence = config["GENOMICS"]["MIN_NB_READS_IN_SEQUENCE"]
min_fraction_of_nb_cpg_in_read = config["GENOMICS"]["MIN_FRACTION_OF_NB_CPG_IN_READ"]
sort_reads_randomly = config["GENOMICS"]["SORT_READS_RANDOMLY"]
nb_cpg_for_padding = config["GENOMICS"]["NB_CPG_FOR_PADDING"]

# Feature prep (Kernel functions)
kernel_fm_nb_values = config["FEATURE_PREP"]["KERNEL_FM_NB_VALUES"]
kernel_fm_bandwidth = config["FEATURE_PREP"]["KERNEL_FM_BANDWIDTH"]
kernel_cov_nb_max = config["FEATURE_PREP"]["KERNEL_COV_NB_MAX"]
kernel_cov_nb_step = config["FEATURE_PREP"]["KERNEL_COV_NB_STEP"]
kernel_cov_bandwidth = config["FEATURE_PREP"]["KERNEL_COV_BANDWIDTH"]
kernel_type = config["FEATURE_PREP"]["KERNEL_TYPE"]

# Define the path to the JSON credentials file
credentials_path = "/appuser/.config/gcloud/application_default_credentials.json"

# Check if the JSON credentials file exists
if os.path.exists(credentials_path):
    # Set the GOOGLE_APPLICATION_CREDENTIALS environment variable to the path of the JSON key file
    os.environ["GOOGLE_APPLICATION_CREDENTIALS"] = credentials_path
    # Assuming 'project' is already defined somewhere in your script
    os.environ["GOOGLE_CLOUD_PROJECT"] = project
    samples_dataset = "samples_250bp"

# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()


# Retrieve Job-defined env vars
BATCH_TASK_INDEX = int(os.getenv("BATCH_TASK_INDEX", 0))
ml_dataset = os.getenv("ML_DATASET")
nb_files_per_task = int(os.getenv("NB_FILES_PER_TASK", 0))
samples_dataset = os.getenv("SAMPLES_DATASET")


def create_schema_fields(variables_dict):
    """Create schema fields from a dictionary of variable names and types."""
    schema_fields = [
        bigquery.SchemaField(name, type_, mode="REQUIRED")
        for name, type_ in variables_dict.items()
    ]
    # logging.info(f"Schema fields created: {schema_fields}")
    return schema_fields


def add_record_field(schema_fields, record_name, fields, mode="REPEATED"):
    """Add a record field to an existing schema."""
    record_field = bigquery.SchemaField(record_name, "RECORD", mode=mode, fields=fields)
    return schema_fields + [record_field]


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


def calculate_cpg_distances(cpgs):
    """
    Calculates the consecutive distances between CpG positions in a list.
    Parameters:
    - cpgs: an array of fm, cov, pos as a dictionary.
    Returns:
    - list of int: A list of distances between consecutive CpG positions.
    """
    # Initialize an empty list to store distances between consecutive CpG positions
    cpg_positions = [int(d["pos"]) for d in cpgs]
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


def generate_feature_arrays(
    row,
    min_nb_reads_in_sequence,
    min_fraction_of_nb_cpg_in_read,
    sort_reads_randomly,
    nb_cpg_for_padding,
):
    reads = row["reads"]
    nb_cpg_found = int(row["nb_cpg_found"])
    # Obtain the reads that have enough CpGs
    reads_w_enough_cpgs = [
        read_info
        for read_info in reads
        if len(read_info["cpg_pos"])
        >= int(np.round(min_fraction_of_nb_cpg_in_read * nb_cpg_found))
    ]
    # Return None if we do not have enough reads
    if len(reads_w_enough_cpgs) < min_nb_reads_in_sequence:
        # print("Returning None")
        return pd.Series([None, None])
    # Ensure the data type of each dictionary is correct
    for read in reads_w_enough_cpgs:
        # Convert 'fm' to float
        read["fm"] = float(read["fm"])
        # Convert 'cpg_pos' to list of integers
        read["cpg_pos"] = [int(pos) for pos in read["cpg_pos"]]
        # Convert 'cpg_meth' to list of integers
        read["cpg_meth"] = [int(meth) for meth in read["cpg_meth"]]
    # Sort reads by FM descending
    reads_sorted_by_fm = sorted(
        reads_w_enough_cpgs, key=lambda d: d["fm"], reverse=True
    )
    # Pick reads at random if variable is set to true and there are more reads than what we need
    if sort_reads_randomly and len(reads_sorted_by_fm) > min_nb_reads_in_sequence:
        # Randomly sample profiles if more are available than the minimum required
        reads_sorted_by_fm = random.sample(reads_sorted_by_fm, min_nb_reads_in_sequence)
    # Find the list of unique CpG positions
    cpg_unique_pos = sorted(
        set(pos for read in reads_sorted_by_fm for pos in read["cpg_pos"])
    )
    # Create a dictionary of all CpGs positions (sorted by position) with their methylation status (0 if no information, 1 for unmethylated CpG, 2 for methylated CpG)
    cpg_dic = {cpg_pos: [] for cpg_pos in cpg_unique_pos}
    meth_count = 0
    for read in reads_sorted_by_fm:
        meth_count += 1
        for cpg_pos, cpg_meth in zip(read["cpg_pos"], read["cpg_meth"]):
            cpg_dic[cpg_pos] += [cpg_meth + 1]
        for unique_cpg in cpg_unique_pos:
            if len(cpg_dic[unique_cpg]) < meth_count:
                cpg_dic[unique_cpg] += [0]
    # Compute directional fractional methylation
    cpg_directional_fm = [
        np.round(
            (
                sum(x[: int(min_nb_reads_in_sequence / 2)])
                - sum(x[int(min_nb_reads_in_sequence / 2) :])
            )
            / min_nb_reads_in_sequence,
            3,
        )
        for x in cpg_dic.values()
    ]
    # Add padding to the cpg_dic
    cpgs_w_padding = [cpg_methyl for cpg_methyl in cpg_dic.values()]
    # First scenario: remove extra CpGs from the trailing end
    while len(cpgs_w_padding) > nb_cpg_for_padding:
        cpgs_w_padding.pop(0)
        if len(cpgs_w_padding) > nb_cpg_for_padding:
            cpgs_w_padding.pop(-1)
    # Second scenario: add CpGs with zeros on each side consecutively
    while len(cpgs_w_padding) < nb_cpg_for_padding:
        cpgs_w_padding = [[0] * min_nb_reads_in_sequence] + cpgs_w_padding
        if len(cpgs_w_padding) < nb_cpg_for_padding:
            cpgs_w_padding += [[0] * min_nb_reads_in_sequence]
    # Return resultats
    return pd.Series([cpg_directional_fm, cpgs_w_padding])


# Define main script
def main():

    logging.info(f"Importing: {nb_files_per_task} files...")
    # logging.info(f"Config file : {config}")

    # Store the JSON file into a dataframe
    df, file_name = create_df_from_json_for_index_file(
        storage_client,
        bucket,
        samples_dataset + "/all_regions/",
        BATCH_TASK_INDEX,
        nb_files_per_task,
    )

    # Uncomment this for testing
    # df_raw = df_raw.head(10)

    logging.info(f"Number of rows in raw dataframe: {len(df)}")

    # logging.info("Create kernel functions")
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

    # logging.info("Calculate consecutive distances between CpGs")
    df["cpg_dist"] = df["cpgs"].apply(calculate_cpg_distances)
    df["cpg_cov"] = df["cpgs"].apply(lambda x: [int(d["cov"]) for d in x])
    df["cpg_fm"] = df["cpgs"].apply(lambda x: [int(d["fm"]) for d in x])
    df["read_fm"] = df["reads"].apply(lambda x: [int(d["fm"]) for d in x])
    df["nb_reads"] = df["read_fm"].apply(len)

    # logging.info("Convert specific arrays into kernel densities")
    for col in dic_kernel.keys():
        logging.info(f"Processing kernel density for {col}")
        if col in df.columns:
            convert_arrays_to_kernel_densities(
                df,
                col,
                dic_kernel[col]["values"],
                dic_kernel[col]["bandwidth"],
            )

    for col in dic_kernel.keys():
        expand_array_elements(df, col + "_kd")

    logging.info(
        "Create two new columns, one with the sequence of directional CpG FM, one with the 2d picture nb CpGs x nb_reads"
    )

    df[
        [
            "cpg_directional_fm",
            "cpgs_w_padding",
        ]
    ] = df.apply(
        lambda x: generate_feature_arrays(
            x,
            min_nb_reads_in_sequence,
            min_fraction_of_nb_cpg_in_read,
            sort_reads_randomly,
            nb_cpg_for_padding,
        ),
        axis=1,
        result_type="expand",
    )

    logging.info("Removing extra columns")
    df = df.drop(
        columns=vars_to_remove + ["cpgs", "reads"], axis=1, errors="ignore"
    ).copy(deep=True)

    logging.info("One-hot encode categoricals variables that are not binary")
    dummies_list = []
    for var in categorical_vars_ohe:
        # logging.info(f"One hot for variable {var}")
        dummies = pd.get_dummies(df[var], prefix=var, dtype=int)
        dummies_list.append(dummies)
    df = pd.concat([df, pd.concat(dummies_list, axis=1)], axis=1)

    for chr_num in range(1, 23):
        if "chr_" + str(chr_num) not in df.columns:
            df["chr_" + str(chr_num)] = 0

    # logging.info("Enforcing data types for integer variables")
    for var in [
        "chr",
        "region_inf",
        "region_sup",
        "clustering_index",
        "region_nb_cpg",
        "nb_reads",
        "nb_cpg_found",
        "asm",
    ]:
        df[var] = df[var].astype(pd.Int64Dtype())

    df[["cpg_directional_fm", "cpgs_w_padding"]] = df[
        ["cpg_directional_fm", "cpgs_w_padding"]
    ].applymap(lambda x: f'"{x}"' if x != "" else x)

    logging.info("Uploading the data to BigQuery")
    upload_dataframe_to_bq(bq_client, df, f"{ml_dataset}.features_wo_hmm")

    logging.info("SCRIPT COMPLETE")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
