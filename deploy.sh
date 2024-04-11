#!/bin/bash

# you need this package to parse the YAML file
# brew install yq

SHORT_SHA="$(git rev-parse --short HEAD)"
echo "SHORT_SHA: ${SHORT_SHA}"

# Submit the build to Google Cloud Build
gcloud builds submit --config=cloudbuild.yaml . --substitutions=SHORT_SHA="${SHORT_SHA}"



# Run locally for faster debugging
# docker run -it \
# -v ~/.config/gcloud/application_default_credentials.json:/root/.config/gcloud/application_default_credentials.json:ro \
# -e GOOGLE_APPLICATION_CREDENTIALS=/root/.config/gcloud/application_default_credentials.json \
# us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:dc47b73 \
# /bin/bash

# # Write this first:
# export GOOGLE_CLOUD_PROJECT=your-project-id
