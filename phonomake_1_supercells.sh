#!/bin/bash

phonopy --qe -d --dim="2 2 2" -c scf.in
#phonopy --qe -d --dim="2 2 2" -c scf.in --tolerance=1e-3

echo ""
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo "  REMEMBER TO COMMENT  !celldm(1)=xxx  FROM  scf.in  "
echo "  AND TO UPDATE nat WITH THE NEW NUMBER OF ATOMS     "
echo "  BEFORE INSERTING THE HEADERS                       "
echo "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
echo ""
