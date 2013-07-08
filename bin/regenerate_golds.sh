#!/bin/bash

GOLDDIR=/raid_hum/piernik_golds

trunk-mcrtest_new() {
   git checkout 7e18cf6dad72004847876079c25251e68fe3f1d
   python setup mcrtest --param problem.par.build2 --compiler gnu47
   pushd runs/mcrtest &> /dev/null
   mpiexec --mca btl ^openib -n 1 piernik -n '&OUTPUT_CONTROL run_id="gld" /'
   mv mcr_gld_0020.h5 ${GOLDDIR}
   popd &> /dev/null
}

trunk-mcrwind_validation() {
   git checkout 7e18cf6dad72004847876079c25251e68fe3f1d
   python setup mcrwind --param problem.par.build --compiler gnu47 --debug
   pushd runs/mcrwind &> /dev/null
   mpiexec --mca btl ^openib -n 4 ./piernik -n '&OUTPUT_CONTROL run_id="gld" use_v2_io= .true. /'
   mv crwind_final_gld_0001.h5 ${GOLDDIR}
   popd &> /dev/null
}

trunk-resist() {
   git checkout 7e18cf6dad72004847876079c25251e68fe3f1d
   python setup tearing --param problem.par.build --compiler gnu47 --debug
   pushd runs/tearing &> /dev/null
   mpiexec --mca btl ^openib -n 1 ./piernik -n '&OUTPUT_CONTROL run_id="gld" use_v2_io=.true. /'
   mv tearing_final_gld_0001.h5 ${GOLDDIR}
   popd &> /dev/null
}

trunk-si() {
   git checkout 1a944b96f840531f51dc638a302da42666841369
   python setup streaming_global --param problem.par.build --compiler gnu47
   pushd runs/streaming_global &> /dev/null
   mpiexec --mca btl ^openib -n 4 ./piernik -n '&OUTPUT_CONTROL run_id="gld" use_v2_io=.true. /'
   mv si_gld_0001.h5 ${GOLDDIR}
   popd &> /dev/null
}

TMPDIR=$(mktemp -d)
pushd ${TMPDIR} &> /dev/null

git clone http://github.com/piernik-dev/piernik
cd piernik

trunk-mcrtest_new
trunk-mcrwind_validation
trunk-resist
trunk-si

popd &> /dev/null
rm -rf ${TMPDIR}
