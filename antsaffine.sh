#!/bin/bash 

NUMPARAMS=$#

if [ $NUMPARAMS -lt 4 ]
then 
echo "USAGE :
antsaffine.sh  fixed.nii  moving.nii  outputvolname.nii  rigidonly(yes/no) numcpu 

FIXED          Fixed image (nifti)
MOVING         Moving image (nifti)
OUTPUTVOL      Output image (nifti .nii)
RIGIDONLY      A flag (yes or no), indicating if the registration is rigid (yes) or affine (no)
NUMCPU         Number of parallel processing cores to be used. Optional."
exit
fi



dim=3
FSLOUTPUTTYPE=NIFTI
START=$(date +%s)
f=$1
m=$2
output=$3
isrigid=$4
numcpu=$5
prefix=`remove_ext $output`
red=`tput setaf 5`
green=`tput setaf 2`
reset=`tput sgr0`

if [ x"$numcpu" == "x" ];then
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=12
else
    ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$numcpu
fi    
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS    

AP=`which antsRegistration`
AP=`readlink -f $AP`
if [ -f "$AP" ];then
    AP=`dirname $AP`
    echo "${green}I found ANTs installation at $AP $reset"
else
    echo "${red}I did not find ANTs in your path. Please install ANTs and add the bin directory to your PATH. $reset"    
    exit 1
fi 

its=100x50x25
percentage=0.1
if [ -f "$prefix".nii ];then
    echo "${red}The following file exists. I will not overwrite."
    echo "$prefix".nii
    echo "If you want to rerun, please delete/move this file first. Exiting.${reset}"
    exit 1
fi   
if [ -f "$prefix".nii.gz ];then
    echo "${red}The following file exists. I will not overwrite."
    echo "$prefix".nii.gz
    echo "If you want to rerun, please delete/move this file first. Exiting.${reset}"
    exit 1
fi  
if [ "$isrigid" == "yes" ];then
    echo "antsRegistration -d $dim -r [ $f, $m ,1]  -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t translation[ 0.1 ] -c [$its,1.e-8,20] -s 4x2x1vox  -f 6x4x2 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t rigid[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -o [ ${prefix} ] -v 0"
    antsRegistration -d $dim -r [ $f, $m ,1]  -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t translation[ 0.1 ] -c [$its,1.e-8,20] -s 4x2x1vox  -f 6x4x2 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t rigid[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -o [ ${prefix} ] -v 0
    
else
    echo "antsRegistration -d $dim -r [ $f, $m ,1]  -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t translation[ 0.1 ] -c [$its,1.e-8,20] -s 4x2x1vox  -f 6x4x2 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t rigid[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t affine[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -o [ ${prefix} ] -v 0"
    antsRegistration -d $dim -r [ $f, $m ,1]  -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t translation[ 0.1 ] -c [$its,1.e-8,20] -s 4x2x1vox  -f 6x4x2 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t rigid[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t affine[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -o [ ${prefix} ] -v 0
    
fi    

antsApplyTransforms -d $dim -i $m -r $f -o ${prefix}.nii -n BSpline -t "$prefix"0GenericAffine.mat -f 0 -v 0
fslmaths ${prefix}.nii -thr 0 ${prefix}.nii

END=$(date +%s)
DIFF=$(( $END - $START ))
((sec=DIFF%60, DIFF/=60, min=DIFF%60, hrs=DIFF/60))
if [ "$isrigid" == "yes" ];then
    echo "${green}ANTS rigid registration took $hrs HRS $min MIN $sec SEC${reset}"
else    
    echo "${green}ANTS affine registration took $hrs HRS $min MIN $sec SEC${reset}"
fi    

