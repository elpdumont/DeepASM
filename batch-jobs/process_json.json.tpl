{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "ML_DATASET": "ML_DATASET_PH",
            "CLOUDASM_DATASET": "CLOUDASM_DATASET_PH",
            "NB_FILES_PER_TASK": "NB_FILES_PER_TASK_PH"
          }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PH",
              "commands": ["python", "/app/process_cloudasm_to_ml_datasets.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 1500,
          "memoryMib": 20000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "2000s"
      },
      "taskCount": TASK_COUNT_PH,
      "parallelism": 60
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": { "machineType": "n1-highmem-96"
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
