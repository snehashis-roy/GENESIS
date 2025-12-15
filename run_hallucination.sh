#!/bin/sh

exe_name=$0
exe_dir=`dirname "$0"`

if [ $# -lt "16" ]; then
echo "----------------------------------------------------------------------------------------------"
echo "Usage:
./run_hallucination.sh atlasvollist subjectvollist isskull NumAtlasPatches
NumSubjectPatches NumNearestNeighbors PatchSize NumCPU OutlierLevel InputModality OutputModality
AtlasStartSlice SubjectStartSlice AtlasWMPeaks SubjectWMPeaks outputvolname MCRROOT
Input and outvolume name should be either NIFTI or XML
atlasvollist          Comma separated atlas volumes, e.g. to synthesize FLAIR from T1 and T2,
                      use spgr.nii,t2.nii,flair.nii
subjectvollist        Comma separated subject volumes, e.g. t1.nii,t2.nii (number of subject
                      volumes = number of atlas volumes - 1)
isskull               Enter yes or no
NumAtlasPatches       Number of atlas patches used to generated atlas patch cloud. A typical
                      number is between 5E4 to 5E5
NumSubjectPatches     Number of subject patches used to generated subject patch cloud.
                      A typical number is between 5E4 to 1E5
NumNearestNeighbors   Number of nearest atlas patches to be searched for each subject patch.
                      A typical number is 15 to 25.
PatchSize             Either 3x3x1 or 3x3x3
NumCPU                For multi-core processors. Default is maximum number of available cpus
Outlierlevel          To remove outliers, either low,medium,high, or none (default)
InputModality         Options are SPGR,MPRAGE,T2 or CT. For multichannel input, use the
                      modality of first channel
OutputModality        Options are SPGR,MPRAGE,T2 or CT
AtlasStartSlice       Atlas slices where the patches are collected from. 0-->Default.
                      Usually choose slices from the middle of the atlas image
SubjectStartSlice     Subject slices where the patches are collected from. 0-->Default.
                      Usually choose slices from the middle of the subject image
AtlasWMPeaks          For images with skull, enter WM peaks in comma separated format,
                      e.g. 199.0,299.0
SubjectWMPeaks        For images with skull, enter WM peaks in comma separated format,
                      e.g. 199.0,299.0
outputvolname         Full path and name of the output volume. Either XML or NIFTI

----------------------------------------------------------------------------------------------"
  exit 1
else
    ATLAS=$1
    SUBJECT=$2
    SKULL=$3
    NUMATLAS=$4
    NUMSUBJECT=$5
    NUMNN=$6
    PSIZE=$7
    NUMCPU=$8
    OUTLIER=$9
    INPUTMODAL=${10}
    OUTPUTMODAL=${11}
    ATLASSLICE=${12}
    SUBJECTSLICE=${13}
    ATLASWMPEAK=${14}
    SUBJECTWMPEAK=${15}
    OUTVOL=${16}

    export MCRROOT=/usr/local/matlab-compiler/v912
    LD_LIBRARY_PATH=.:${MCRROOT}/runtime/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/bin/glnxa64 ;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/os/glnxa64;
    LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${MCRROOT}/sys/opengl/lib/glnxa64;
    export LD_LIBRARY_PATH;
  

    echo "${exe_dir}"/run_hallucination $ATLAS $SUBJECT $SKULL $NUMATLAS $NUMSUBJECT $NUMNN $PSIZE $NUMCPU $OUTLIER $INPUTMODAL $OUTPUTMODAL $ATLASSLICE $SUBJECTSLICE $ATLASWMPEAK $SUBJECTWMPEAK $OUTVOL
    "${exe_dir}"/run_hallucination  $ATLAS $SUBJECT $SKULL $NUMATLAS $NUMSUBJECT $NUMNN $PSIZE $NUMCPU $OUTLIER $INPUTMODAL $OUTPUTMODAL $ATLASSLICE $SUBJECTSLICE $ATLASWMPEAK $SUBJECTWMPEAK $OUTVOL

fi
