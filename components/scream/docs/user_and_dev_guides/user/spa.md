# Generating Prescribed Aerosol Climatology (in EAMf90):

By default the SCREAM compsets are run with prescribed aerosols, which is counter to E3SM which is run with prognostic aerosols.  To run with prescribed aerosols requires climatological aerosol output derived from a prognostic aerosol run.  

The functionality to output the relevant data to generate the files needed to drive a prescribed aerosol run is now a feature in SCREAM F90 code.  **NOTE:** it is NOT necessary to generate these files each time a prescribed aerosol run is performed.  The prescribed aerosol input files only need to be generated when substantial modifications are made to the model (i.e. major physics changes or changes to the vertical grid).  

Should the prescribed aerosol files need to be regenerated the user should follow these steps:

1. Set up a global simulation with **prognostic** aerosols.
2. In the namelist (user_nl_cam) be sure to set presc_aero_data = .true.  
3. Run the simulation for six years minimum.  There is a minimum of five years needed to generate the prescribed aerosol climatology and we toss out the first year to allow the aerosols to spin up.  
4. Note that the required aerosol data will be placed in the monthly \*.cam.h0\*.nc files.
5. When the run has completed, additional post processing is needed to get the output data in the structure required for the prescribed aerosol code.  To achieve this, run the script found on the repo scream-docs at design-docs/prescribed_aerosol/make_presc_aero.csh .  Follow the instructions in the header to set the appropriate settings.  

# Running Prescribed Aerosol Simulation (in EAMf90):

When the updated prescribed aerosol file has been generated, the user will need to update their namelist (user_cam_nl) to reflect the updated location and name of the newly generated prescribed aerosol file.  Specifically, the namelist settings that need to be changed are:

1. aerodep_flx_datapath and prescribed_aero_datapath : These will be set to the location of the newly generated prescribed aerosol file (both will be set to the same value).
2. aerodep_flx_file and prescribed_aero_file : These will be set to the file name of the newly generated aerosol file (both will be set to the same value).   

Then the user should be able to run with any supported FSCREAM compset.  
