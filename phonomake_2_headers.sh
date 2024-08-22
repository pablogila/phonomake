#!/bin/bash

NAT_BY_HAND=0
HEADER_FILE='scf.in'
KEYWORD='!BEGIN_COORDINATES'

# Check if the supercell contains '&CONTROL' and abort if it does
supercelltest=$(ls supercell-*.in | head -n 1)
if grep -q '&CONTROL' "$supercelltest"; then
  echo "ERROR: Found '&CONTROL' in $supercelltest, so seems like you already did this!"
  return
fi

# Extract the content from scf.in up to the keyword indicating the start of the old coordinates
awk "/$KEYWORD/ {exit} {print}" "$HEADER_FILE" > temp_header.txt
echo "Extracted header from $HEADER_FILE"

# Comment all 'celldm(1)' and 'nat' lines except those that are already commented
sed -i '/^\s*[^!]\s*celldm(1)\s*=/s/\(celldm(1)\s*=\)/!\1/' temp_header.txt
sed -i '/^\s*[^!]\s*nat\s*=/s/\(nat\s*=\)/!\1/' temp_header.txt
echo "Commented old celldm(1) and nat values"

# Update nat with the number of atoms in the supercell, or manually if not found
nat=$(grep -oP 'nat\s*=\s*\K\d+' "$supercelltest")
if [ -n "$nat" ]; then
    sed -i "0,/!nat =/!b; /!nat =/a\  nat = $nat" temp_header.txt
    echo "Updated  nat = $nat"
else
    if [ "$NAT_BY_HAND" -ne 0 ]; then
        sed -i "0,/!nat =/!b; /!nat =/a\  nat = $NAT_BY_HAND" temp_header.txt
        echo "Updated  nat = $NAT_BY_HAND  manually, from NAT_BY_HAND inside this script"
    else
        echo "ERROR: No 'nat' value found in $supercell."
        echo "       Set it manually inside this script"
        echo "       with  NAT_BY_HAND=nat"
        rm temp_header.txt
        return
    fi
fi

# Loop through all supercell-*.in files and prepend the new header
for file in supercell-*.in
do
    temp_file=$(mktemp)
    cat temp_header.txt "$file" > "$temp_file" && mv "$temp_file" "$file"
    echo "Added header to $file"
done

# Remove the temporary file
rm temp_header.txt
