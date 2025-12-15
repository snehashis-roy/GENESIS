#!/bin/sh
exe_name=$0
exe_dir=`dirname "$0"`
if [ $# -lt "3" ]; then
  echo "---------------------------------------------------------------------"
  echo "Usage:"
  echo "./run_remove_background_noise.sh INPUT QUANTILE OUTPUT"
  echo "All volumes must be either NIFTI or XML"
  echo "INPUT       : An MRI 3D volume, with skull and with background noise"
  echo "QUANTILE    : Percent of histogram to consider for initial quantile.
              A reasonable value is 20-40. (default 30)"
  echo "OUTPUT      : Final noise removed volume"
  echo "---------------------------------------------------------------------"
  exit 1
else
    INPUT=$1 
    Q=$2
    OUTPUT=$3
    MCRROOT=/usr/local/matlab-compiler/v912
    LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
    export LD_LIBRARY_PATH;

    echo "${exe_dir}"/remove_background_noise $INPUT $Q $OUTPUT 
    "${exe_dir}"/remove_background_noise $INPUT $Q $OUTPUT 
fi
exit

