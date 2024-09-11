# PhonoMake

Easy and automated workflow for Phonopy calculations using Quantum ESPRESSO.


## Workflow

Copy the `/pm/` directory into your calculation folder, as well as the `scf.slurm` file.  

We start from a relaxed, geometry optimized structure. To use the *PhonoMake* scripts, we need to have two keywords on the `relax.in` input file:
- `!END_HEADER`, which goes right after `K_POINTS`, and before `ATOMIC_SPECIES`, `CELL_PARAMETERS` and `ATOMIC_POSITIONS`
- `!BEGIN_COORDINATES`, which goes right after `ATOMIC_SPECIES`, before `CELL_PARAMETERS` and `ATOMIC_POSITIONS`

```relax.in
&CONTROL
&SYSTEM
...
K_POINTS
!END_HEADER
ATOMIC_SPECIES
!BEGIN_COORDINATES
CELL_PARAMETERS
ATOMIC_POSITIONS
```

Once relaxed, we create a `scf.in` file by running:
```shell
source pm/0_scf.sh
```

Check that the resulting `scf.in` file is properly configured, eg. by running the `cat scf.in` command.  

Then, create the 2x2 supercells with
```shell
source pm/1_supercells.sh
```

Or if you want a different supercell, edit the shell command in the `pm/1_supercells.sh` file.  

Next, add the headers (taken from the `scf.in` file) with
```shell
source pm/2_headers.sh
```

Cat a random supercell to check that everything is ok.  

By now, you should have copied the provided `scf.slurm` file inside your folder.
Modify it to your needs, but do not change the `JOB_NAME`, `INPUT_FILE` nor `OUTPUT_FILE` keywords.  

To submit the calculations, run
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


## What are these scripts doing to my beloved files?

### pm/0_scf.sh

The `pm/0_scf.sh` script performs the following tasks.  

To create the `scf.in` file, it first copies the input parameters from `relax.in`, up to the `!BEGIN_COORDINATES` keyword, right after `ATOMIC_SPECIES` and before `CELL_PARAMETERS` and `ATOMIC_POSITIONS`.  

Then, it copies the `CELL_PARAMETERS` and `ATOMIC_POSITIONS` sections from the `relax.out` file.  

Phonopy [requires](https://phonopy.github.io/phonopy/qe.html) a `&SYSTEM celldm(1)` flag to specify the lattice parameter in bohr units. Notice that the previous lattice parameter `A` introduced by [cif2cell](https://github.com/torbjornbjorkman/cif2cell) is in Armstrongs. We can find the proper celldm(1) value in the CELL_PARAMETERS output (alat= 16.7...). The CELL_PARAMETERS line is then left as `CELL_PARAMETERS alat`, such as follows:

```scf.in
&SYSTEM 
	!A = 8.86370  ! Old lattice parameter used to relax
	celldm(1) = 16.74996545  ! Add this line with the proper alat value
	...

!CELL_PARAMETERS (alat= 16.74996545)  ! Leave only 'alat', as follows:
CELL_PARAMETERS alat
	...
```

### pm/1_supercells.sh

It just creates the supercells, by running the phonopy command:
```shell
phonopy --qe -d --dim="2 2 2" -c scf.in
```

### pm/2_headers.sh

It extracts the header from `scf.in` up to the `!END_HEADER` keyword, and adds it to the new supercells. It also updates the number of atoms, and comments the old lattice parameter which is no longer needed.

### pm/3_slurms.sh
It creates the slurm files for each supercell from the provided `scf.slurm` template, and sbatch'es the corresponding jobs.
You can modify `scf.slurm` to your needs, just try not to change the `JOB_NAME`, `INPUT_FILE` and `OUTPUT_FILE` keywords.

### pm/4_forces.sh

It creates the `FORCE_SETS` file from the QE output of the supercell calculations,
and uses it to create the `phonopy.yaml` file, by running
```shell
phonopy -f supercell-*.out
phonopy --include-all
```

### pm/fix_scancel.sh
Scancels all the jobs with a `slurm-JOB_ID` file in the folder.

### pm/fix_unfinished.sh
Resubmits all unfinished supercells. It asks for a new memory value, in case the previous one was not enough.

### pm/fix_yaml.sh
In case Phonopy does not recognise some of your atomic species and changes it by a random one, you can use this to fix the final `phonopy.yaml` file.
However, this is already [fixed](https://github.com/phonopy/phonopy/issues/412) if your custom atomic species follow a **symbol+number** format, so luckily you will never have to use this script.
