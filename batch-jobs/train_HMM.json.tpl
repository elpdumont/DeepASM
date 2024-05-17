{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "ML_DATASET": "ML_DATASET_PH",
            "MODEL_PATH": "MODEL_PATH_PH",
            "SHORT_SHA": "SHORT_SHA_PH"
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
          "cpuMilli": 96000,
          "memoryMib": 360000
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
        "policy": { "machineType": "n1-standard-96"
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
