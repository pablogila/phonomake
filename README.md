# PhonoMake

Workflow for Phonopy calculations using Quantum ESPRESSO.

## Workflow

Copy the `/pm/` directory into your calculation folder, as well as the `scf.slurm` file.  

We start from a relaxed, geometry optimized structure. Once relaxed we have to create a `scf.in` file, copying the final coordinates from the relaxed output file.  

Phonopy [requires](https://phonopy.github.io/phonopy/qe.html) a `&SYSTEM celldm(1)` flag to specify the lattice parameter in bohr units. Notice that the previous lattice parameter `A` introduced by [cif2cell](https://github.com/torbjornbjorkman/cif2cell) is in Armstrongs. We can find the proper celldm(1) value in the CELL_PARAMETERS output (alat= 16.7...). The CELL_PARAMETERS line must be then left as `CELL_PARAMETERS alat`, such as follows:

```scf.in
&SYSTEM 
	! A = 8.86370  ! Old lattice parameter used to relax
	celldm(1) = 16.74996545  ! Add this line with the proper alat value
	...

! CELL_PARAMETERS (alat= 16.74996545)  ! Leave only 'alat', as follows:
CELL_PARAMETERS alat
	...
```

To use the *PhonoMake* scripts, you also need to  separate the coordinates below `K_POINTS` by a `!BEGIN_COORDINATES` keyword. That way, the headers can be introduced automatically.  

Check that the `scf.slurm` file is properly configured, and you are ready to start the calculations.  

First, create the 2x2 supercells with
```shell
source pm/1_supercells.sh
```

If you want a different supercell, edit the shell command in the `pm/1_supercells.sh` file.  

Then, add the headers with
```shell
source pm/2_headers.sh
```

Cat a random supercell to check that everything is ok. To submit the calculations, use
```shell
source pm/3_slurms.sh
```

Finally, once the calculations are done, calculate the forces with
```shell
source pm/4_forces.sh
```

In addition, there are some scripts to fix common errors in the calculations:
- `pm/fix_scancel.sh`
- `pm/fix_unfinished.sh`
- `pm/fix_yaml.sh`
