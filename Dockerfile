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

# Set the default entry point to Python
#ENTRYPOINT ["python"]

CMD ["python", "./process_json.py"]

# To test locally with Docker
# docker pull us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest
# docker run --rm -it us-east1-docker.pkg.dev/hmh-em-deepasm/docker-repo/python:latest python3 --version