#!/bin/bash
# Ben Fasoli
set -e

echo "Executing setup script..."
./setup

echo "Sourcing R dependencies..."
Rscript r/dependencies.r

echo "Setup successful."
