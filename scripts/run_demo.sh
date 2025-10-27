#!/usr/bin/env bash
set -e
echo "Running iDNDS demo (dN/dS calculation)..."

if [ -d ".venv" ]; then
  source .venv/bin/activate
fi

python src/main.py

echo "Demo finished successfully."
echo "Output file: data/expected_output/demo_result.txt"
