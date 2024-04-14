{
  "taskGroups": [
    {
      "taskSpec": {
        "environment": {
          "variables": {
            "ML_DATASET_ID": "ML_DATASET_ID_PLACEHOLDER"
            }
        },
        "runnables": [
          {
            "container": {
              "imageUri": "PYTHON_IMAGE_PLACEHOLDER:IMAGE_TAG_PLACEHOLDER",
              "commands": ["python", "/app/run_hmm.py"]
            }
          }
        ],
        "maxRetryCount": 0,
        "maxRunDuration": "20000s"
      },
      "taskCount": 1,
      "parallelism": 1
    }
  ],
  "allocationPolicy": {
    "instances": [
      {
        "policy": { "machineType": "e2-highmem-8"
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