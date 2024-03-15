# File management
import json
import logging
import os
import random
import sys

import dask
import dask.dataframe as dd

# Python packages for data, stats
import numpy as np
import pandas as pd
import yaml
from google.cloud import bigquery, storage
from sklearn.neighbors import KernelDensity

from gcp_utils import create_df_from_json_for_index_file

# Initialize random seed (for selecting the reads used in the matrix)
random_seed = 42
random.seed(random_seed)


# Initialize the Google Cloud Storage client
storage_client = storage.Client()
bq_client = bigquery.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# Disable warning log for dask
dask.config.set({"dataframe.query-planning-warning": False})

# Retrieve Job-defined env vars
TASK_INDEX = int(os.getenv("CLOUD_RUN_TASK_INDEX", 0))
TASK_ATTEMPT = os.getenv("CLOUD_RUN_TASK_ATTEMPT", 0)

with open("config.yaml", "r") as file:
    config = yaml.safe_load(file)

# Accessing GCP configuration
bucket_name = config["GCP"]["BUCKET_NAME"]
raw_data_bucket_folder = config["GCP"]["RAW_DATA_BUCKET_FOLDER"]
bq_ml_dataset_name = config["GCP"]["BQ_ML_DATASET_NAME"]

# Variables to handle
vars_to_remove = config["VARS_TO_REMOVE"]
dic_vars_to_keep = config["VARS_TO_KEEP"]
vars_to_keep = [list(item.keys())[0] for item in dic_vars_to_keep]

# vars_to_normalize = config["VARS_TO_NORMALIZE"]
categorical_vars = config["CATEGORICAL_VARS"]
categorical_vars_ohe = config["CAT_VARS_OHE"]

# Genomic variables
genomic_length = config["GENOMICS"]["GENOMIC_LENGTH"]
min_cpg_cov = config["GENOMICS"]["MIN_CPG_COV"]

# Feature prep (Kernel functions)
kernel_fm_nb_values = config["FEATURE_PREP"]["KERNEL_FM_NB_VALUES"]
kernel_fm_bandwidth = config["FEATURE_PREP"]["KERNEL_FM_BANDWIDTH"]
kernel_cov_nb_max = config["FEATURE_PREP"]["KERNEL_COV_NB_MAX"]
kernel_cov_nb_step = config["FEATURE_PREP"]["KERNEL_COV_NB_STEP"]
kernel_cov_bandwidth = config["FEATURE_PREP"]["KERNEL_COV_BANDWIDTH"]
kernel_type = config["FEATURE_PREP"]["KERNEL_TYPE"]


def compute_counts_and_percentages(column):
    """
    Computes the counts and percentages of zeros and ones in a column.

    Args:
    - column (pd.Series): A pandas Series containing the data.

    Returns:
    - counts (dict): A dictionary with counts of zeros and ones.
    - percentages (dict): A dictionary with percentages of zeros and ones.
    """
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
            lambda x: x[i] if i < len(x) else None
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
    import numpy as np

    # Initialize the list to hold dictionaries for each genomic position
    result = [{"pos": k + 1, "cpg": 0, "cpg_fm": None} for k in range(genomic_length)]

    # Update dictionaries in the result list for positions with CpGs
    for pos, fm in zip(row["cpg_pos"], row["cpg_fm"]):
        # Adjust for 0-based indexing and update only the necessary keys
        index = pos - 1  # Adjusting for 0-based indexing
        result[index]["cpg"] = 1
        result[index]["cpg_fm"] = np.round(fm, 3)

    return result

    # Initialize the list to hold dictionaries for each genomic position
    result = [{f"{k+1}": 0, "cpg_fm": 0.0} for k in range(genomic_length)]

    # Update dictionaries in the result list for positions with CpGs
    for pos, fm in zip(row["cpg_pos"], row["cpg_fm"]):
        result[pos] = {f"{pos}": 1, "cpg_fm": fm}

    return result


def generate_methyl_profile_for_read(read_dictionary):
    pos_array = np.array(read_dictionary["pos_array"], dtype=int)
    meth_array = np.array(read_dictionary["meth_array"], dtype=int)

    # Initialize the result array with zeros
    result = np.zeros((genomic_length, 2), dtype=int)

    # Set the first column to 1 at positions where CpG sites are present
    result[pos_array - 1, 0] = 1  # Assuming pos_array is 1-indexed

    # Update the methylation status
    result[pos_array - 1, 1] = meth_array

    return {
        "read_id": read_dictionary["read_id"],
        "read_fm": read_dictionary["read_fm"],
        "nb_cpg": read_dictionary["nb_cpg"],
        "cpg_methylation_profile": result,
    }


def generate_cpg_methylation_matrix(row_df, min_coverage):
    """
    Generates a 3D CpG methylation matrix for a set of genomic reads
    within a specific genomic window.
    The matrix dimensions are length_genomic_interval x min_coverage x 2, where:
    - The first dimension (length_genomic_interval) represents the length of the genomic window.
    - The second dimension (min_coverage) corresponds to the minimum number of reads required.
    - The third dimension contains a tuple [a, b] for each genomic position in a read:
        'a' indicates the presence (1) or absence (0) of a CpG site,
        'b' indicates the methylation status of the site if present.

    Reads are ordered by fractional methylation (FM) in descending order within the matrix.
    If the minimum coverage is not met, an empty list is returned.

    Parameters:
    - row_df (pd.Series): A row from a DataFrame, representing data for a specific genomic window.
    - genomic_interval (int): The length of the genomic window to be considered.
    - min_coverage (int): The minimum number of reads required to generate the matrix.

    Returns:
    - methylation_matrix (np.array or list):
      A 3D numpy array representing the CpG methylation matrix if the minimum coverage is met;
      otherwise, an empty list. For instance, if the genomic window length is 250 and the min coverage is 20,
      the shape of the matrix will be (250, 20, 2).
    """

    genomic_array = row_df["genomic_picture"]

    # Generate the 1D CpG methylation profile for each read
    profiles = [
        generate_methyl_profile_for_read(read_dict) for read_dict in genomic_array
    ]

    # Filter out empty results
    profiles = [profile for profile in profiles if profile]

    if profiles and len(profiles) >= min_coverage:

        # Randomly sample profiles if more are available than the minimum required
        sampled_profiles = (
            random.sample(profiles, min_coverage)
            if len(profiles) > min_coverage
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

        # Transpose the matrix to match the required dimensions
        methylation_matrix = np.transpose(methylation_matrix, (1, 0, 2))
    else:
        methylation_matrix = []

    return np.array(methylation_matrix)


# Define main script
def main():

    # Store the JSON file into a dataframe
    df_raw, file_name = create_df_from_json_for_index_file(
        bucket_name, raw_data_bucket_folder, TASK_INDEX
    )
    logging.info(f"File name: {file_name}")
    logging.info(f"Number of rows in raw dataframe: {len(df_raw)}")

    logging.info("Create kernel functions")
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

    logging.info("Remove unwanted chromosomes")
    df_filtered = filter_chromosomes(df_raw)
    logging.info(f"Number of rows of df_filtered: {len(df_filtered)}")
    logging.info(compute_counts_and_percentages(df_filtered["asm_snp"]))

    logging.info("Calculate consecutive distances between CpGs")
    df_filtered["cpg_dist"] = df_filtered["cpg_pos"].apply(calculate_cpg_distances)

    logging.info("Convert specific arrays into kernel densities")
    for col in dic_kernel.keys():
        if col in df_filtered.columns:
            convert_arrays_to_kernel_densities(
                df_filtered,
                col,
                dic_kernel[col]["values"],
                dic_kernel[col]["bandwidth"],
            )

    logging.info(
        "Create a column for each kernel density estimate (right now they are stored in arrays)"
    )
    for col in dic_kernel.keys():
        expand_array_elements(df_filtered, col + "_kd")

    logging.info(
        "Create a methylation sequence of nucleotide over the genomic length. Each nucleotide comes with 2 infos: presence of CpG, presence of methylation"
    )
    df_filtered["sequence_cpg_fm"] = df_filtered.apply(
        lambda x: generate_nucleotide_sequence_of_cpg_fm(x, genomic_length), axis=1
    )

    logging.info(
        "Create a methylation matrix of nucleotide over the genomic length. Each nucleotide is represented by a vector of the depth. Each element represents the presence of a CpG and its methylation status"
    )
    ddf = dd.from_pandas(df_filtered, npartitions=10)
    result = ddf.apply(
        lambda x: generate_cpg_methylation_matrix(x, min_cpg_cov),
        axis=1,
        meta=("x", "object"),
    )
    df_filtered["methylation_matrix"] = result.compute()

    initial_row_count = len(df_filtered)
    # Remove rows where the specified column contains an empty list
    df_filtered = df_filtered[
        df_filtered["methylation_matrix"].apply(lambda x: len(x) > 0)
    ]

    # Count the number of rows after removing rows with empty lists
    final_row_count = len(df_filtered)

    # Calculate the number of rows removed
    rows_removed = initial_row_count - final_row_count

    logging.info(
        f"There are {rows_removed} rows without methylation matrix data from the dataset out of {initial_row_count} rows"
    )
    percentage_removed = (rows_removed / initial_row_count) * 100
    logging.info(f"{percentage_removed:.2f}% of the rows were removed")

    # Store the different datasets into a hash table.
    dic_data = {}

    dic_data["clean"] = df_filtered.drop(
        columns=vars_to_remove, axis=1, errors="ignore"
    )

    logging.info("One-hot encode categoricals variables that are not binary")
    dic_data["clean"] = pd.get_dummies(
        dic_data["clean"], columns=categorical_vars_ohe, dtype=int
    )

    logging.info("Create a list of categorical variables")
    categorical_vars = [
        col
        for col in dic_data["clean"].columns
        if any(col.startswith(var) for var in categorical_vars)
    ]

    logging.info(f"All categorical variables: {categorical_vars}")

    logging.info("Creating the dataset with a methylation sequence")
    dic_data["sequence_cpg_fm"] = dic_data["clean"][
        vars_to_keep + ["sequence_cpg_fm"]
    ].copy(deep=True)

    logging.info("Creating the  dataset with the methylation matrix")
    dic_data["methylation_matrix"] = dic_data["clean"][
        vars_to_keep + ["methylation_matrix"]
    ].copy(deep=True)

    logging.info("Creating the tabular dataset")
    dic_data["tabular"] = (
        dic_data["clean"]
        .drop(columns=["sequence_cpg_fm", "methylation_matrix"], axis=1)
        .copy(deep=True)
    )

    logging.info("Uploading the 3 datasets to BigQuery")
    schema_fields = []

    for name, type_ in dic_vars_to_keep.items():
        # Here we assume all fields are required, adjust if that's not the case
        field = bigquery.SchemaField(name, type_, mode="REQUIRED")
        schema_fields.append(field)

    record_fields_sequence_cpg_fm = [
        bigquery.SchemaField("pos", "INTEGER", mode="NULLABLE"),
        bigquery.SchemaField("cpg", "INTEGER", mode="NULLABLE"),
        bigquery.SchemaField("cpg_fm", "FLOAT", mode="NULLABLE"),
    ]
    record_field_sequence_cpg_fm = bigquery.SchemaField(
        "sequence_cpg_fm",
        "RECORD",
        mode="REPEATED",
        fields=record_fields_sequence_cpg_fm,
    )
    schema_sequence_cpg_fm = schema_fields.append(record_field_sequence_cpg_fm)

    logging.info(
        "Uploading the dataframe with the sequence of CpG fractional methylations"
    )
    job_config = bigquery.LoadJobConfig(schema=schema_sequence_cpg_fm)
    bq_client.load_table_from_dataframe(
        dic_data["sequence_cpg_fm"], "ml.test7", job_config=job_config
    )

    # for data_type in ["tabular", "sequence_cpg_fm", "methylation_matrix"]:
    #     export_df_to_gcs_json(
    #         dic_data[data_type],
    #         data_type,
    #         BUCKET_NAME,
    #         OUTPUT_DATA_FOLDER_PATH,
    #         file_name,
    #     )


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = (
            f"Task #{TASK_INDEX}, " + f"Attempt #{TASK_ATTEMPT} failed: {str(err)}"
        )

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
