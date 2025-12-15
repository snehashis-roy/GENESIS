#!/bin/sh
exe_name=$0
exe_dir=`dirname "$0"`
if [ $# -lt "2" ]; then
  echo "---------------------------------------------------------------------"
  echo "Usage:
  ./find_WM_peak.sh  input_image  modality [brainmask]
 input_image       Either stripped or unstripped T1/T2/PD/FL image. If the
                   image is unstripped, a brain mask is required for 
                   correct estimation of white matter peak
 modality          T1/T2/PD/FL
 brainmask         (Optional) Not required if the input image is stripped. For images
                   with skull, it is required.
  
  "
  exit 1
fi


MCRROOT=/usr/local/matlab-compiler/v912
LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
export LD_LIBRARY_PATH;

args=
while [ $# -gt 0 ]; do
  token=$1
  args="${args} ${token}" 
  shift
done
"${exe_dir}"/find_WM_peak $args

