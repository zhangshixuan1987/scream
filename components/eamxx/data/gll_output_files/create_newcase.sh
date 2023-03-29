#!/bin/bash

# this file must be *sourced*
[[ "${BASH_SOURCE[0]}" == "${0}" ]] && echo "this script must be sourced" && exit

source ~/mappy-module-list-cime.sh

ttasks=16

scream_top_dir=/ascldap/users/tccleve/E3SM/SCREAM/scream       #or wherever you cloned the scream repo
COMPSET=F2010-SCREAMv1      #or whatever compset you want
RES=ne4pg2_ne4pg2         #or whatever resolution you want
CASE_NAME=Test-GLL-Output-${COMPSET}-${RES}-${ttasks}ranks  #name for your simulation.
casepath=${CIME_SCRATCH_DIR}/${CASE_NAME}
COMPILER=gnu9             #which compiler to use (can be omitted on some machines)
MACH=mappy

cd ${scream_top_dir}/cime/scripts/
./create_newcase -mach ${MACH} --case ${casepath} --compset ${COMPSET} --res ${RES} --compiler ${COMPILER}  --handle-preexisting-dirs r

cd ${casepath}

./xmlchange NTASKS=${ttasks} #change how many MPI tasks to use
./xmlchange DEBUG=TRUE #debug build.
./xmlchange STOP_OPTION=nsteps #how long to run for
./xmlchange STOP_N=1
./xmlchange PIO_NETCDF_FORMAT="64bit_data"
./xmlchange SCREAM_CMAKE_OPTIONS="`./xmlquery -value SCREAM_CMAKE_OPTIONS | sed 's/SCREAM_NUM_VERTICAL_LEV [0-9][0-9]*/SCREAM_NUM_VERTICAL_LEV 128/'`"
./xmlchange HIST_N=1; ./xmlchange HIST_OPTION=nsteps ;
./xmlchange REST_N=999999; ./xmlchange REST_OPTION=nyears ;

./case.setup

#has to be done afgter setup
./atmchange disable_diagnostics=true
./atmchange output_yaml_files=${scream_top_dir}/components/eamxx/data/gll_output_files/output_with_gll.yaml

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

ffol=tmp-run${id}
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
