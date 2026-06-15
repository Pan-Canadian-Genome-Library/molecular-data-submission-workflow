FROM node:22-slim

LABEL org.opencontainers.image.source="https://github.com/Pan-Canadian-Genome-Library/molecular-data-submission-workflow"

RUN apt-get update && apt-get install -y --no-install-recommends \
    procps \
    && rm -rf /var/lib/apt/lists/*

# Copy and install dependencies
COPY package*.json ./
RUN npm install --omit=dev
