# Use an official Python runtime as a parent image
FROM python:3.10.12-slim

# Create a non-root user
RUN adduser --disabled-password --gecos '' appuser

# Set the working directory in the container
WORKDIR /app

# First, copy only the requirements.txt file to leverage Docker cache
COPY requirements.txt .

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the current directory contents into the container at /app
COPY app/gcp-functions-python /app
COPY app/prepare-samples/python /app

# Copy the config files
COPY config.yaml /app

# Switch to the non-root user
USER appuser

CMD ["python"]

# To test locally with Docker
# docker pull us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:f59e8b7
# docker run --rm -it --platform linux/amd64 us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:f59e8b7

# Then, in the image, type:
# import os
# print(os.listdir())
