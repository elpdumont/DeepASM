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
              "commands": ["python", "find_nb_states_for_HMM.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 47000,
          "memoryMib": 180000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "500000s"
      },
      "taskCount": TASK_COUNT_PH,
      "parallelism": 2
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
