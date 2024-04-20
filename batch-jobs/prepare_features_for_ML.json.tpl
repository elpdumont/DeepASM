{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "ML_DATASET": "ML_DATASET_PH",
            "SAMPLES_DATASET": "SAMPLES_DATASET_PH",
            "NB_FILES_PER_TASK": "NB_FILES_PER_TASK_PH"
          }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PH",
              "commands": ["python", "prepare_features_for_ML.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 3500,
          "memoryMib": 30000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "5000s"
      },
      "taskCount": TASK_COUNT_PH
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
