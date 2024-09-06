#!/bin/bash

phonopy --qe -d --dim="2 2 2" -c scf.in
#phonopy --qe -d --dim="2 2 2" -c scf.in --tolerance=1e-3
