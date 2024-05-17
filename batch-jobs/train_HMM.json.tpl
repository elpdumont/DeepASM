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
              "commands": ["python", "train_HMM.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 120000,
          "memoryMib": 860000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "500000s"
      },
      "taskCount": 1
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": { "machineType": "n2-highmem-128"
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
