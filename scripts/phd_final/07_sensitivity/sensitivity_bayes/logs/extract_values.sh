#!/usr/bin/env bash

counter=1

for file in *.log; do

    echo "Processing $file"

    # --- extract signal process ---
    signal=$(grep -m1 "signal process:" "$file" | awk -F': ' '{print $2}')

    # AAA = first token before underscore
    AAA=$(echo "$signal" | cut -d'_' -f1)

    # --- extract RADON_TAG ---
    radon=$(grep -m1 "RADON_TAG:" "$file" | awk -F': ' '{print $2}')

    # --- build output filename ---
    outfile="${AAA}_tag_${radon}_${counter}.csv"

    echo "Writing $outfile"

    # --- extract sensitivities ---
    # print line AFTER "time to fit"
    awk '
        /time to fit/ {
            getline
            print
        }
    ' "$file" > "$outfile"

    ((counter++))

done
