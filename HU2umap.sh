#!/bin/sh
exe_name=$0
exe_dir=`dirname "$0"`
if [ $# -lt "3" ]; then
  echo "---------------------------------------------------------------------"
  echo "Usage:"
  echo "./run_HU2umap.sh INPUTFILE OUTPUTFILE KVP"
  echo "All volumes must be either NIFTI or XML"
  echo "INPUTFILE         A CT 3D volume"
  echo "OUTPUTFILE        Mu-map from CT"
  echo "KVP               Default kVp 120"
  echo "---------------------------------------------------------------------"
  exit 1
else
  INPUT=$1 
  OUTPUT=$2
  KVP=$3
  MCRROOT=/usr/local/matlab-compiler/v912
  LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
  MCRJRE=${MCRROOT}/sys/java/jre/glnxa64/jre/lib/amd64 ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/native_threads ; 
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/server ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE}/client ;
  LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRJRE} ;  
  XAPPLRESDIR=${MCRROOT}/X11/app-defaults ;
  export LD_LIBRARY_PATH;
  export XAPPLRESDIR;
  echo "${exe_dir}"/HU2umap $INPUT $OUTPUT  $KVP
  "${exe_dir}"/HU2umap $INPUT $OUTPUT $KVP
fi
exit

