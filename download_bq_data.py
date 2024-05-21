import json
import os

from google.cloud import bigquery

# Initialize BigQuery client
client = bigquery.Client()

# Define the datasets
datasets = ['hmh-em-deepasm.ml_250bp_db7e6e4', 'hmh-em-deepasm.ml_250bp_70efde8']

# Create a directory to store the JSON files
output_dir = 'bigquery_json_output'
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

def download_table_to_json(dataset_id, table_id, output_dir):
    table_ref = client.dataset(dataset_id).table(table_id)
    table = client.get_table(table_ref)

    rows = client.list_rows(table)
    rows_list = [dict(row) for row in rows]
    
    file_path = os.path.join(output_dir, f"{dataset_id}_{table_id}.json")
    with open(file_path, 'w') as json_file:
        json.dump(rows_list, json_file, indent=2)

    print(f"Downloaded {table_id} to {file_path}")

for dataset in datasets:
    dataset_id = dataset.split('.')[-1]
    tables = client.list_tables(dataset)
    for table in tables:
        download_table_to_json(dataset_id, table.table_id, output_dir)

print("All tables downloaded successfully.")