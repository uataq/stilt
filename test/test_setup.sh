#!/bin/bash
# Ben Fasoli
set -e

# Execute setup script to install stilt
./setup

# Ensure R package dependencies are installed
Rscript r/dependencies.r

echo "Setup successful."
