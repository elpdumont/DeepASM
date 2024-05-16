{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "ML_DATASET": "ML_DATASET_PH",
            "MODEL_PATH": "MODEL_PATH_PH",
            "SHORT_SHA": "SHORT_SHA_PH",
            "ML_MODE": "ML_MODE_PH"
            }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PH",
              "commands": ["python", "perform_random_search_tree.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 10000,
          "memoryMib": 80000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "60000s"
      },
      "taskCount": 2,
      "parallelism": 2
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": { "machineType": "n1-highmem-32"
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
