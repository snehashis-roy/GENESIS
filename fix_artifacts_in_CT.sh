#!/bin/bash

exe_name=$0
exe_dir=`dirname "$0"`

if [ $# -lt "3" ]; then
echo "----------------------------------------------------------------------------"
echo "Usage:
./run_fix_artifacts_in_CT.sh CT BRAINMASK MODIFIED_CT
CT               CT image, generated using run_hallucination.sh
BRAINMASK        Brain mask, usually generated using robex/spectre
MODIFIED_CT      Modified CT image, where errors inside the brainmask are fixed
----------------------------------------------------------------------------"
exit 1
else
    CT=$1  
    MASK=$2
    MODCT=$3
    export MCRROOT=/usr/local/matlab-compiler/v912
    LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
    export LD_LIBRARY_PATH;
          
    echo "${exe_dir}"/fix_artifacts_in_CT $CT $MASK $MODCT
   "${exe_dir}"/fix_artifacts_in_CT $CT $MASK $MODCT
   
fi
