#!/bin/bash

# Phonopy has an issue (#412) in which it replaces 'unknown' atomic species, such as `H2`, to random ones
# If the pseudos are properly set, it should not affect the previous Quantum ESPRESSO SCF calculation
# This script should be a quick fix for the final phonopy.yaml
# AbINS reads `H`, and knows that are deuteriums because of the mass, etc.

FILE=phonopy.yaml

echo "Let's fix the atomic species that Phonopy messed up in the $FILE file!"
read -p "Symbol for the old atomic species to replace:  " OLD_ATOMS
read -p "Symbol for the new atomic species:  " NEW_ATOMS

KEYWORD_OLD="  - symbol: $OLD_ATOMS "
KEYWORD_NEW="  - symbol: $NEW_ATOMS "

tmpfile=$(mktemp)
sed "s|$KEYWORD_OLD|$KEYWORD_NEW|g" "$FILE" > "$tmpfile"
if cmp -s "$FILE" "$tmpfile"; then
    echo "No match found, no changes were made."
    rm "$tmpfile"
else
    mv "$tmpfile" "$FILE"
    echo "Replaced keys in $FILE"
    echo "$KEYWORD_OLD >>> $KEYWORD_NEW"
fi
