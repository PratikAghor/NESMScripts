* ```NESM_2019_2020, NESM_2019_2020_5km, SM3_2019_2020``` contain scripts to run and plot post-processed data for the NESM and SEAMOUNT (SM3) runs, respectively.
* ```figs``` directory contains data (wherever possible) and scripts for plotting figures in the manuscript. 
* ```aghor_extras``` directory has other packages required to plot figures, see README.md there for details.

### Steps to configure simulations and postprocess output

#### Download bulk and boundary data from HYCOM and NAVGEM: `blk_bry_data`
* In directory `blk_bry_data/HYCOM`, use `gethycom_varname.sbatch` to download the individual HYCOM files,
* Then use `blk_bry_data/merged/merge_hycom_data.sbatch` to consolidate individual HYCOM ncfiles into one file per variable
* For NAVGEM data, follow instructions in the README file at `blk_bry_data/NAVGEM`.
* Tidal forcing and GOT99.2 can be added from CROCO's website to `blk_bry_data/TPXO7` and `blk_bry_data/GOT99.2`

#### Location of Gridfiles
* Gridfile for NESM (1 km) run: `NESM_2019_2020/data/edit_grid/NESM_grd.nc`. 
* Gridfile for NESM (5 km) run: `NESM_2019_2020_5km/data/edit_grid/NESM_grd_5km.nc`. 
* Gridfile for SM3 (1 km) run: `SM3_2019_2020/data/edit_grid/SM3_grd.nc`.
* Note that these gridfiles must be copied to the `forcing` directory before running the simulations.
* Create a directory to store forcings `mkdir forcing`.
* Copy gridfiles from `data/edit_grid/` to the `forcing` directory.

#### Generating forcing files and initial conditions: `forcing_prep`
* Each run contains a `forcing_prep` directory, where all the forcings (boundary, tidal, etc.) are interpolated to be consistent with the corresponding grid.
* One must adjust `dir_path`, `blk_bry_path`according to the system in the `make_BLK_navgem.m`, `make_BRY_hycom.m`, 'make_FRC_tides.m` and `make_INI_hycom.m`.
* Also adjust paths in `start_paths.m`.
* Use `job_make_BLK.sbatch` to generate the '../forcing/*_blk*_.nc` (bulk) file. 
* Use `job_make_BRY.sbatch` to generate the '../forcing/*_OGCM*_.nc` and '../forcing/*_bry*_.nc` (boundary) files. 
* Use `job_make_FRC_tides.sbatch` to generate the '../forcing/*_frc_tides_*.nc` (tideal forcing) file.
* Use `job_make_INI.sbatch` to generate the '../forcing/*_ini*_.nc` (initial conditions) file.
* At the end of this step, the `forcing` directory should contain six ncfiles: `_grd.nc`, `*_blk_*.nc`, `*_OGCM_*.nc`, `*_bry_*.nc`, `*_frc_tides_*.nc` and `*_ini_*.nc`. 

#### Compiling and running the simulations: `compile`
* Before this step, generate an output directory where output will be saved: `mkdir data/output`.
* Now switch to the compile directory `cd compile`.
* Change paths accordingly in the `compile/croco.in` file for forcings to the appropriate `forcing` directory and also update paths for saving the output to the appropriate `data/output` directory. 
* Update `param.h` if you need to configure the run with different number of processors, etc.
* `cppdefs.h` contains details for numerical schemes, best not to update it to reproduce the same results. 
* `./jobcomp > comp.log` to compile the run, it should generate a `coroco` executable. 
* Finally, use `sbatch job_*.sbatch` for running the simulation. 
* Make sure that the data is being saved in the `data/output` directory correctly. 

#### Postprocessing: `postprocessing_scripts` and `figs`
* Once the `data/output` is generated, use matlab scripts in `postprocessing_scripts` to postprocess the data.
* Adjust paths everywhere in mostly `start_paths.m` according to your system. 
* `postprocessing_scripts/aghor_Akt`: Contains scripts related to vertical diffusivity ($k_v$) calculations. 
* `postprocessing_scripts/aghor_eke`: Contains scripts related to eddy kinetic energy (EKE) calculations.
* `postprocessing_scripts/aghor_kmke`: Contains scripts related to KmKe calculations.
* `postprocessing_scripts/aghor_ow`: Contains scripts related to Okubo-Weiss calculations.
* `postprocessing_scripts/aghor_pv`: Contains scripts related to Ertel potential vorticity ($q$) calculations.
* `postprocessing_scripts/aghor_vebf`: Contains scripts related to VEBF (vertical eddy buoyancy flux) calculations.
* `postprocessing_scripts/aghor_vot`: Contains scripts related to relative vorticity ($\zeta$) calculations.
* `postprocessing_scripts/aghor_w`: Contains scripts related to vertical velocity (w) calculations.  
* Finally, `figs` directory contains scripts and data to generate figures in the manuscript.     