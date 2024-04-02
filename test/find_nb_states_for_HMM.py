import json
import logging
import os
import sys

import numpy as np
import pandas as pd
import yaml
from gcp.utils import upload_blob
from google.cloud import bigquery, storage
from hmmlearn.hmm import GaussianHMM
from sklearn.utils import check_random_state


def main():
    logging.info("Starting script")


# Start script
if __name__ == "__main__":
    try:
        main()
    except Exception as err:
        message = f"Task #{BATCH_TASK_INDEX} failed: {str(err)}"

        print(json.dumps({"message": message, "severity": "ERROR"}))
        sys.exit(1)  # Retry Job Task by exiting the process
