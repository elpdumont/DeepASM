{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "SAMPLES_DATASET": "SAMPLES_DATASET_PH"
            }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PH",
              "commands": ["python", "calculate_wilcoxon_for_regions.py"]
            }
          }
        ],
        "computeResource": {
          "cpuMilli": 5000,
          "memoryMib": 20000
        },
        "maxRetryCount": 0,
        "maxRunDuration": "2000s"
      },
      "taskCount": NB_SAMPLES_PH,
      "parallelism": NB_SAMPLES_PH
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
