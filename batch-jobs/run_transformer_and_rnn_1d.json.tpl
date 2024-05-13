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
              "commands": ["python", "run_transformer_and_rnn_1d.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 10000,
          "memoryMib": 30000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "12000000s"
      },
      "taskCount": 2,
      "parallelism": 2
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "installGpuDrivers": true,
        "policy": { "machineType": "g2-standard-24"
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
