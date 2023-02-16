#!/bin/bash

# this file must be *sourced*
[[ "${BASH_SOURCE[0]}" == "${0}" ]] && echo "this script must be sourced" && exit

source ~/summit-modules-cime.sh

####switch to ranks
### nodes 4 8,16,32,43,65,129,300, (600 does not work)
nnodes=32
pernode=6
ttasks=$(expr $pernode \* $nnodes)

scream_top_dir=~/E3SM-Project/SCREAM/scream       #or wherever you cloned the scream repo
COMPSET=F2010-SCREAMv1-X6T      #or whatever compset you want
RES=ne30_ne30         #or whatever resolution you want
CASE_NAME=${COMPSET}-${RES}-default-settings-${nnodes}nodes  #name for your simulation.
casepath=${CIME_SCRATCH_DIR}/${CASE_NAME}
COMPILER=gnugpu             #which compiler to use (can be omitted on some machines)
MACH=summit

cd ${scream_top_dir}/cime/scripts/
./create_newcase -mach ${MACH} --case ${casepath} --compset ${COMPSET} --res ${RES} --compiler ${COMPILER}  --handle-preexisting-dirs r

cd ${casepath}

./xmlchange NTASKS=${ttasks} #change how many MPI tasks to use
./xmlchange MAX_TASKS_PER_NODE=6
./xmlchange MAX_MPITASKS_PER_NODE=6
#./xmlchange DEBUG=TRUE #debug build.
./xmlchange DEBUG=FALSE #optimized build.
./xmlchange JOB_QUEUE=batch #debug or batch
./xmlchange JOB_WALLCLOCK_TIME=2:00
./xmlchange STOP_OPTION=nmonths #how long to run for
./xmlchange STOP_N=14
./xmlchange PIO_NETCDF_FORMAT="64bit_data"
./xmlchange SCREAM_CMAKE_OPTIONS="`./xmlquery -value SCREAM_CMAKE_OPTIONS | sed 's/SCREAM_NUM_VERTICAL_LEV [0-9][0-9]*/SCREAM_NUM_VERTICAL_LEV 128/'`"
./xmlchange HIST_N=1; ./xmlchange HIST_OPTION=nmonths ;
./xmlchange REST_N=999999; ./xmlchange REST_OPTION=nyears ;

./case.setup

#has to be done afgter setup
./atmchange disable_diagnostics=true
./atmchange pgrad_correction=0
./atmchange hv_ref_profiles=0
./atmchange output_yaml_files=${scream_top_dir}/components/eamxx/data/rough_topo_tests_output.yaml

if true; then

./case.build
./case.submit

#do it properly
#./xmlquery EXEROOT
#	EXEROOT: ${casepath}/bld
bbld=${casepath}/bld

#e3sm
#id=`grep -o case.run:\[0-9\]\* CaseStatus | sed -e 's/.*://' | tail -1`
#scream
id=`grep -o case.run:\[0-9\]\* env_run.xml  | sed -e 's/.*://' | tail -1`

echo $id

ffol=run${id}
mkdir $ffol
cd $ffol
cp ../namelist_scream.xml namelist_scream."${id}".xml
#cp ../user_nl_eam usercam."${id}"
cp ${bbld}/GIT_*.* .

gitstat=gitstat."${id}"
cd ${scream_top_dir}
echo " running stats on clone ${scream_top_dir}"
echo " status is ------------- " >> ${gitstat}
git status >> ${gitstat}
echo " branch is ------------- " >> ${gitstat}
git branch >> ${gitstat}
echo " diffs are ------------- " >> ${gitstat}
git diff >> ${gitstat}
echo " last 10 commits are ------------- " >> ${gitstat}
git log --first-parent  --pretty=oneline  HEAD~10..HEAD  >> ${gitstat}

cd ${casepath}/${ffol}
mv ${scream_top_dir}/${gitstat} .
cd ${casepath}


fi
