#!/bin/bash
set -e

echo "Installing SNP Analyzer dependencies..."

# Upgrade pip
pip install --upgrade pip setuptools wheel

# Install packages
pip install PyQt5
pip install pysam
pip install pandas
pip install matplotlib  
pip install seaborn
pip install biopython

echo "✓ All dependencies installed successfully!"
echo "Running SNP Analyzer..."

python desktop_app.py
