{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "SAMPLE_LIST": "SAMPLE_LIST_PH",
            "SAMPLES_DATASET": "SAMPLES_DATASET_PH",
            "CLOUDASM_DATASET": "CLOUDASM_DATASET_PH",
            "CLOUDASM_TABLES": "CLOUDASM_TABLES_PH"
            }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "BASH_IMAGE_PH",
              "commands": ["/bin/bash", "overlap_cloudasm_data_with_standard_regions.sh"]
            }
          }
        ],
        "maxRetryCount": 0,
        "maxRunDuration": "20000s"
      },
      "taskCount": nb_tasks_for_cloudasm_tables_ph,
      "parallelism": 16
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": { "machineType": "e2-standard-16"
      }}
    ]
  },
  "labels": {
    "env": "deployment"
  },
  "logsPolicy": {
    "destination": "CLOUD_LOGGING"
  }
}
