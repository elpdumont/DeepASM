# # Kernel functions
import json
import logging

# File and variable management
import os
import random
import sys

import dask
import dask.dataframe as dd

# Python packages for data, stats
import numpy as np
import pandas as pd
from google.cloud import storage
from numpy import asarray, exp
from sklearn.neighbors import KernelDensity

# Initialize random seed
random_seed = 42
random.seed(random_seed)
np.random.seed(random_seed)


# Initialize the Google Cloud Storage client
storage_client = storage.Client()

# Initialize logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s")

# Disable warning log for dask
dask.config.set({"dataframe.query-planning-warning": False})


# Retrieve Job-defined env vars
TASK_INDEX = int(os.getenv("CLOUD_RUN_TASK_INDEX", 0))
TASK_ATTEMPT = os.getenv("CLOUD_RUN_TASK_ATTEMPT", 0)

# Retrieve User-defined env vars
BUCKET_NAME = os.getenv("BUCKET_NAME")
INPUT_DATA_FOLDER_PATH = os.getenv("INPUT_DATA_FOLDER_PATH")
OUTPUT_DATA_FOLDER_PATH = os.getenv("OUTPUT_DATA_FOLDER_PATH")
GENOMIC_LENGTH = int(os.getenv("GENOMIC_LENGTH"))
MIN_CPG_COV = int(os.getenv("MIN_CPG_COV"))
KERNEL_FM_NB_VALUES = int(os.getenv("KERNEL_FM_NB_VALUES"))
KERNEL_FM_BANDWIDTH = float(os.getenv("KERNEL_FM_BANDWIDTH"))
KERNEL_COV_NB_MAX = int(os.getenv("KERNEL_COV_NB_MAX"))
KERNEL_COV_NB_STEP = int(os.getenv("KERNEL_COV_NB_STEP"))
KERNEL_COV_BANDWIDTH = int(os.getenv("KERNEL_COV_BANDWIDTH"))

VARS_TO_REMOVE = [
    "region_inf",
    "region_sup",
    "read_fm",
    "cpg_pos",
    "cpg_cov",
    "cpg_dist",
    "cpg_fm",
    "genomic_picture",
    "global_cpg_fm",
    "tot_nb_reads",
    "tot_nb_cpg",
    "cpg_cov_kd",
    "read_fm_kd",
    "cpg_fm_kd",
    "cpg_dist_kd",
]

VARS_TO_NORMALIZE = [
    "region_nb_cpg",
    "nb_cpg_found",
    "nb_reads",
    "encode_ChiP_V2",
    "tf_motifs",
    "cpg_cov_kd",
    "read_fm_kd",
    "cpg_fm_kd",
    "cpg_dist_kd",
    "std_read_fm",
    "mean_read_fm",
    "std_cpg_fm",
    "mean_cpg_fm",
    "std_cpg_cov",
    "mean_cpg_cov",
    "std_cpg_dist",
    "mean_cpg_dist",
]

CATEGORICAL_VARS = ["chr", "sample_category", "dnase"]

CATEGORICAL_VARS_TO_ONE_HOT_ENCODE = ["chr"]


def export_df_to_gcs_json(df, data_type, bucket_name, bucket_path, file_name):
    # Convert DataFrame to JSON format
    json_data = df.to_json(orient="records", lines=True)

    file_name = data_type + "_" + file_name
    logging.info(f"Saving JSON data as {file_name}")
    with open(file_name, "w") as file:
        file.write(json_data)

    # Cloud filename includes the bucket_path
    cloud_filename = os.path.join(bucket_path, file_name)

    # Initialize a Google Cloud Storage client
    storage_client = storage.Client()

    # Get the bucket
    bucket = storage_client.bucket(bucket_name)

    # Create a blob (GCS object) and upload the file
    blob = bucket.blob(cloud_filename)
    blob.upload_from_filename(file_name)

    logging.info(f"File {file_name} uploaded to {cloud_filename}.")


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


def estimate_kernel_density(x, values_for_kernel, bandwidth, kernel="gaussian"):
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
    logging.info(f"Calculating the standard deviation and mean for {column_name}")
    df[std_name] = df[column_name].apply(
        lambda x: np.round(np.std(x, ddof=1), 4)
    )  # ddof=1 for sample standard deviation
    df[mean_name] = df[column_name].apply(lambda x: np.round(np.mean(x), 4))

    # Kernel density estimates
    kernel_name = column_name + "_kd"
    logging.info(f"Calculating the probability distribution for {column_name}")
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


def generate_methylation_sequence(row, genomic_length=GENOMIC_LENGTH):
    """
    Creates arrays of positions and fractional methylation for CpGs within a specified genomic region.
    Positions with CpGs are marked with 1, others with 0. The fractional methylation is set according to the CpG positions,
    with 0 for positions without CpGs.

    Parameters:
    - row pd.Series: Must contain 'cpg_pos', and 'cpg_fm' keys/columns.

    Returns:
    - np.array: An array of shape (n, 2) where 'n' is the length of the genomic region, with the first column indicating
                the presence of CpGs and the second column indicating fractional methylation at those positions.
    """

    # Initialize arrays with zeros for positions and fractional methylation
    genomic_positions = np.zeros(genomic_length, dtype=int)
    genomic_fm = np.zeros(genomic_length)

    # Convert CpG positions to zero-based index relative to the region's start
    # cpg_indices = np.array(row['cpg_pos'].astype(int))

    # Mark positions with CpGs and assign fractional methylation
    genomic_positions[row["cpg_pos"]] = 1
    genomic_fm[row["cpg_pos"]] = row["cpg_fm"]

    return np.stack([genomic_positions, genomic_fm], axis=1)


def generate_methyl_profile_for_read(read_dictionary):
    pos_array = np.array(read_dictionary["pos_array"], dtype=int)
    meth_array = np.array(read_dictionary["meth_array"], dtype=int)

    # Initialize the result array with zeros
    result = np.zeros((GENOMIC_LENGTH, 2), dtype=int)

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


def process_file_at_index(bucket_name, folder_path, task_index):
    # Initialize the GCP Storage client
    storage_client = storage.Client()

    # Define the prefix to search within a specific folder, ensuring it ends with '/'
    prefix = folder_path if folder_path.endswith("/") else f"{folder_path}/"

    # List all blobs in the specified folder
    blobs = storage_client.list_blobs(bucket_name, prefix=prefix)

    # Filter blobs that match the 'raw-*.json' pattern and sort them
    filtered_blobs = sorted(
        (
            blob
            for blob in blobs
            if blob.name.startswith(prefix + "raw-") and blob.name.endswith(".json")
        ),
        key=lambda blob: blob.name,
    )

    # Ensure task_index is within the range of available files
    if 0 <= task_index < len(filtered_blobs):
        # Get the Blob object at the specified index
        selected_blob = filtered_blobs[task_index]
    else:
        logging.info(f"Task index {task_index} is out of range.")
        sys.exit(1)

    # Log the processing info
    logging.info(
        f"Processing the bucket {bucket_name} with the folder path {folder_path} and the file name {selected_blob.name}"
    )

    logging.info("Download the file as bytes and decode it to a string")
    file_contents = selected_blob.download_as_bytes().decode("utf-8")

    logging.info("Process each line as a separate JSON object")
    processed_data = []
    for line in file_contents.splitlines():
        data = json.loads(line)  # Parse each line as JSON
        processed_data.append(data)

    return pd.DataFrame(processed_data), os.path.basename(selected_blob.name)


# Define main script
def main(bucket_name, folder_path, max_digits=12):

    # Store the JSON file into a dataframe
    df_raw, file_name = process_file_at_index(bucket_name, folder_path, TASK_INDEX)
    logging.info(f"File name: {file_name}")
    logging.info(f"Number of rows in raw dataframe: {len(df_raw)}")

    logging.info("Create kernel functions")
    values_for_kernel_cov = generate_kernel_values(
        0, KERNEL_COV_NB_MAX, KERNEL_COV_NB_STEP
    )
    bandwidth_cov = KERNEL_COV_BANDWIDTH

    # Kernel function for fractional methylation ('read_fm', 'cpg_fm')
    values_for_kernel_fm = generate_kernel_values(0, 1, step=1 / KERNEL_FM_NB_VALUES)
    bandwidth_fm = KERNEL_FM_BANDWIDTH

    dic_kernel = {
        "read_fm": {"values": values_for_kernel_fm, "bandwidth": bandwidth_fm},
        "cpg_fm": {"values": values_for_kernel_fm, "bandwidth": bandwidth_fm},
        "cpg_cov": {"values": values_for_kernel_cov, "bandwidth": bandwidth_cov},
        "cpg_dist": {"values": values_for_kernel_cov, "bandwidth": bandwidth_cov},
    }

    logging.info("Remove unwanted chromosomes")
    df_filtered = filter_chromosomes(df_raw)
    logging.info(f"Number of rows of df_filtered: {len(df_filtered)}")
    logging.info(compute_counts_and_percentages(df_filtered["asm_snp"]))

    logging.info("Calculate consecutive distances between CpGs")
    df_filtered["cpg_dist"] = df_filtered["cpg_pos"].apply(calculate_cpg_distances)

    logging.info("Example of an array of cpg distances")
    logging.info(df_filtered.iloc[0]["cpg_dist"])

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
    df_filtered["methylation_sequence"] = df_filtered.apply(
        generate_methylation_sequence, axis=1
    )

    logging.info(
        "Create a methylation matrix of nucleotide over the genomic length. Each nucleotide is represented by a vector of the depth. Each element represents the presence of a CpG and its methylation status"
    )
    ddf = dd.from_pandas(df_filtered, npartitions=10)
    result = ddf.apply(
        lambda x: generate_cpg_methylation_matrix(x, MIN_CPG_COV),
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
        columns=VARS_TO_REMOVE, axis=1, errors="ignore"
    )

    logging.info("One-hot encode categoricals variables that are not binary")
    dic_data["clean"] = pd.get_dummies(
        dic_data["clean"], columns=CATEGORICAL_VARS_TO_ONE_HOT_ENCODE
    )

    logging.info("Create a list of categorical variables")
    categorical_vars = [
        col
        for col in dic_data["clean"].columns
        if any(col.startswith(var) for var in CATEGORICAL_VARS)
    ]

    logging.info(f"All categorical variables: {categorical_vars}")

    logging.info("Creating the dataset with a methylation sequence")
    dic_data["methylation_sequence"] = dic_data["clean"][
        ["sample", "asm_snp", "methylation_sequence"]
    ].copy(deep=True)

    logging.info("Creating the  dataset with the methylation matrix")
    dic_data["methylation_matrix"] = dic_data["clean"][
        ["sample", "asm_snp", "methylation_matrix"]
    ].copy(deep=True)

    logging.info("Creating the tabular dataset")
    dic_data["tabular"] = (
        dic_data["clean"]
        .drop(columns=["methylation_sequence", "methylation_matrix"], axis=1)
        .copy(deep=True)
    )

    logging.info("Uploading the 3 datasets to the bucket")
    for data_type in ["tabular", "methylation_sequence", "methylation_matrix"]:
        export_df_to_gcs_json(
            dic_data[data_type],
            data_type,
            BUCKET_NAME,
            OUTPUT_DATA_FOLDER_PATH,
            file_name,
        )


# Start script
if __name__ == "__main__":
    try:
        main(BUCKET_NAME, INPUT_DATA_FOLDER_PATH)
    except Exception as err:
        message = (
            f"Task #{TASK_INDEX}, " + f"Attempt #{TASK_ATTEMPT} failed: {str(err)}"
        )

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
