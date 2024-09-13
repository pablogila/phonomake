#!/bin/bash

# This script will create the scf.in file for the phonopy calculations.
# It will first read the relax.in file up to the !BEGIN_COORDINATES keyword.
# Then it will append the final coordinates from the relax.out file.
# Finally, it will update the proper lattice parameters and calculation type.

RELAX_IN=relax.in
RELAX_OUT=relax.out
SCF_IN=scf.in
KEYWORD='!BEGIN_COORDINATES'
KEYWORD_HEADER='!END_HEADER'
JOB_DONE="JOB DONE."
TIME_EXCEEDED="Maximum CPU time exceeded"
STEPS_EXCEEDED="The maximum number of steps has been reached."

echo "Creating $SCF_IN file..."

# Check if the relax files exist
if [ ! -f "$RELAX_IN" ]; then
    echo "ERROR: No $RELAX_IN file found!"
    return
fi
if [ ! -f "$RELAX_OUT" ]; then
    echo "ERROR: No $RELAX_OUT file found!"
    return
fi

# Check if relax.out is finished
tail=$(tail -n 200 "$RELAX_OUT")
if [[ $tail != *"$JOB_DONE"* ]]; then
    echo "ERROR: $RELAX_OUT did not finish properly"
    return
elif [[ $tail == *"$TIME_EXCEEDED"* ]]; then
    echo "ERROR: $RELAX_OUT stopped due to time limit"
    return
elif [[ $tail == *"$STEPS_EXCEEDED"* ]]; then
    echo "ERROR: $RELAX_OUT stopped due to max steps limit"
    return
fi

# Check if the keyword is found in the relax file
if ! grep -q "$KEYWORD" "$RELAX_IN"; then
    echo "ERROR: '$KEYWORD' keyword not found in $RELAX_IN"
    echo "       Write the keyword after the ATOMIC_SPECIES, just before"
    echo "       the CELL_PARAMETERS and ATOMIC_POSITIONS"
    return
fi

# Check if the header keyword is found in the relax file, otherwise display a warning
if ! grep -q "$KEYWORD_HEADER" "$RELAX_IN"; then
    echo "WARNING: '$KEYWORD_HEADER' keyword not found in $RELAX_IN"
    echo "         If you plan on using pm/2_headers.sh,"
    echo "         make sure to add it after the K_POINTS, just before"
    echo "         the ATOMIC_SPECIES, CELL_PARAMETERS and ATOMIC_POSITIONS"
fi



# Copy the header from relax.in, up to the !BEGIN_COORDINATES keyword
awk "/$KEYWORD/ {exit} {print}" "$RELAX_IN" > $SCF_IN
echo "Copied input parameters from $RELAX_IN"

# Append final coordinates from relax.out
sed -n '/Begin final coordinates/,/End final coordinates/p' $RELAX_OUT | sed -n '/CELL_PARAMETERS/,/End final coordinates/p' | sed '/End final coordinates/d' >> $SCF_IN
echo "Appended final coordinates from $RELAX_OUT"

# Extract alat lattice parameter from SCF_IN
# It is written in the format CELL_PARAMETERS (alat= 24.43623749)
alat=$(grep -oP 'alat=\s*\K\d+\.\d+' $SCF_IN)
echo "alat = $alat"
echo "Updating input parameters..."

# Comment calculation= line
sed -i '/^\s*[^!]\s*calculation\s*=/s/\(calculation\s*=\)/!\1/' $SCF_IN
# Write calculation='scf' under !calculation= line
sed -i "/^\s*!calculation\s*=/a\  calculation='scf'" $SCF_IN

# Comment the 'A=' line
sed -i '/^\s*[^!]\s*A\s*=/s/\(A\s*=\)/!\1/' $SCF_IN
# Write celldm(1) = $alat under !A= line
sed -i "/^\s*!A\s*=/a\  celldm(1) = $alat" $SCF_IN

# Comment the 'CELL_PARAMETERS' line
sed -i '/^\s*CELL_PARAMETERS/s/^/!/' $SCF_IN
# Add a line with $CELL_PARAMETERS_OUT under the one starting with !CELL_PARAMETERS
sed -i "/^\s*!CELL_PARAMETERS/a\CELL_PARAMETERS alat" $SCF_IN

echo "Done!"
