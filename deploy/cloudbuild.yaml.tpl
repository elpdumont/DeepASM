substitutions:
  _REGION: REGION_PH
  _PROJECT_ID: PROJECT_PH
  _ARTIFACT_REGISTRY_REPO: ARTIFACT_REGISTRY_REPO_PH
  _IMAGE_TAG: "IMAGE_TAG_PH"

steps:
  # Docker Build
  - name: gcr.io/cloud-builders/docker
    args:
      [
        build,
        -f,
        docker/python/Dockerfile,
        -t,
        "${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/python:${_IMAGE_TAG}",
        .,
      ]
    dir: .

  # Docker push to Google Artifact Registry
  - name: gcr.io/cloud-builders/docker
    args:
      [
        push,
        "${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/python:${_IMAGE_TAG}",
      ]

  - name: gcr.io/cloud-builders/docker
    args:
      [
        build,
        -f,
        docker/bash/Dockerfile,
        -t,
        "${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/bash:${_IMAGE_TAG}",
        .,
      ]
    dir: .

  # Docker push to Google Artifact Registry
  - name: gcr.io/cloud-builders/docker
    args:
      [
        push,
        "${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/bash:${_IMAGE_TAG}",
      ]

# Store images in Google Artifact Registry
images:
  - ${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/python:${_IMAGE_TAG}
  - ${_REGION}-docker.pkg.dev/${_PROJECT_ID}/${_ARTIFACT_REGISTRY_REPO}/bash:${_IMAGE_TAG}
