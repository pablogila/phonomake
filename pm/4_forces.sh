#!/bin/bash

# creates FORCE_SETS file from QE output of supercell calculations
phonopy -f supercell-*.out
#phonopy --tolerance 1e-3 --include-all -f supercell-*.out

# Reads the new FORCE_SETS file and creates a phonopy.yaml file
phonopy --include-all
