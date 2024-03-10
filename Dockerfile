# Use an official Python runtime as a parent image
FROM python:3.12

# Set the working directory in the container
WORKDIR /app

# First, copy only the requirements.txt file to leverage Docker cache
COPY requirements.txt .

# Install any needed packages specified in requirements.txt
RUN pip install --no-cache-dir -r requirements.txt

# Copy the current directory contents into the container at /app
COPY src/ .

# Make port 8080 available to the world outside this container
EXPOSE 8080

# Define environment variable
ENV PORT 8080

# Set the default entry point to Python
#ENTRYPOINT ["python"]

# Run process_json.py when the container launches, without CMD specifying script arguments
#CMD ["python3"]

# docker run --rm -it us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/process_json:v12 python3 --version