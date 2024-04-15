FROM google/cloud-sdk:471.0.0-slim

# Update apt cache and install yq
RUN apt-get update && apt-get install -y wget && \
    wget https://github.com/mikefarah/yq/releases/download/v4.27.5/yq_linux_amd64 -O /usr/local/bin/yq && \
    chmod +x /usr/local/bin/yq

# Create a non-root user
RUN adduser --disabled-password --gecos '' appuser

# Set the working directory in the container
WORKDIR /app

# Copy the bash files to be used.
COPY /app/prepare-samples/bash /app

# Make all Bash scripts executable
RUN chmod +x /app/*.sh

# Copy the config files
COPY config.yaml /app

# Switch to the non-root user
USER appuser

CMD ["bash"]
