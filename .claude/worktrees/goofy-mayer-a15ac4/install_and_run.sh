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

# Resolve Qt platform plugin path from whichever PyQt5 the active Python sees.
# Without this macOS throws "Could not find Qt platform plugin cocoa".
QT_PLUGIN_DIR=$(python -c "
import os, PyQt5
base = os.path.dirname(PyQt5.__file__)
for candidate in [
    os.path.join(base, 'Qt5', 'plugins', 'platforms'),
    os.path.join(base, 'Qt',  'plugins', 'platforms'),
]:
    if os.path.isdir(candidate):
        print(candidate)
        break
" 2>/dev/null)

if [ -n "$QT_PLUGIN_DIR" ]; then
    export QT_QPA_PLATFORM_PLUGIN_PATH="$QT_PLUGIN_DIR"
fi

python desktop_app.py
