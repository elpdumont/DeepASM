{
  "taskGroups": [
    {
      "taskSpec": {
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PH",
              "commands": ["python", "/app/prepare-samples/python/calculate_wilcoxon_for_regions.py"]
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
