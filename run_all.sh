#!/bin/bash

echo "Running Zipf analysis..."

for file in input/*.txt
do
    echo "-----------------------------------"
    echo "Processing: $file"

    ./exe/Zipf.exe "$file"
done

echo "-----------------------------------"
echo "All files processed."