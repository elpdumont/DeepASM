FROM google/cloud-sdk:471.0.0-slim

# Create a non-root user
RUN adduser --disabled-password --gecos '' appuser

# Set the working directory in the container
WORKDIR /app

# Copy the bash files to be used.
COPY ../../app/prepare-samples/bash /app

# Copy the config files
COPY ../../config.yaml /app

# Switch to the non-root user
USER appuser

CMD ["bash"]
