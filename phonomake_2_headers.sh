#!/bin/bash

HEADER_FILE='scf.in'
KEYWORD='!BEGIN_COORDINATES'

# Extract the content from scf.in up to the keyword indicating the start of the coordinates
awk "/$KEYWORD/ {exit} {print}" "$HEADER_FILE" > temp_header.txt
echo "Extracted header from $HEADER_FILE"
# Loop through all supercell-*.in files and prepend the extracted content
for file in supercell-*.in
do
    temp_file=$(mktemp)
    cat temp_header.txt "$file" > "$temp_file" && mv "$temp_file" "$file"
    echo "Added header to $file"
done
# Remove the temporary file
rm temp_header.txt
