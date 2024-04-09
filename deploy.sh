#!/bin/bash

# you need this package to parse the YAML file
# brew install yq

SHORT_SHA="$(git rev-parse --short HEAD)"
echo "SHORT_SHA: ${SHORT_SHA}"

# Submit the build to Google Cloud Build
gcloud builds submit --config=cloudbuild.yaml . --substitutions=SHORT_SHA="${SHORT_SHA}"