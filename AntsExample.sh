#!/bin/bash
if [ $# -lt "4" ]; then
echo "==========================================================================
Usage:
$0 fixed.nii moving.nii mysetting registeredvolume.nii

fixed.nii       Fixed image
moving.nii      Moving image
mysetting       Either forproduction (slowest), fast, or fastfortesting(fastest)


========================================================================="
exit 1
fi


red=`tput setaf 5`
green=`tput setaf 2`
reset=`tput sgr0`

dim=3 # image dimensionality


export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS
f=$1 ; m=$2    # fixed and moving image file names
mysetting=$3
outvol=$4


FSLOUTPUTTYPE=NIFTI
prefix=`remove_ext ${outvol}`
reg=antsRegistration           # path to antsRegistration

if [[ $mysetting == "fastfortesting" ]] ; then
  its=100x50x25
  percentage=0.1
  syn="100x1x0,0.0001,4"
elif   [[ $mysetting == "forproduction" ]] ; then
  its=1000x1000x1000
  percentage=0.3
  syn="100x100x50,0.00001,5"
elif [[ $mysetting == "fast" ]] ; then
   its=100x100x100
  percentage=0.3
  syn="20x20x10,0.00001,5"
else
    echo "${red}ERROR: setting must be either forproduction (slowest), fast, or fastfortesting(fastest). $reset"
    exit 1
fi
START=$(date +%s)
echo "$reg -d $dim -r [ $f, $m ,1]  -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t translation[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 6x4x2 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t rigid[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1 -m mattes[  $f, $m , 1 , 32, regular, $percentage ] -t affine[ 0.1 ] -c [$its,1.e-8,20]  -s 4x2x1vox  -f 3x2x1 -l 1   -m mattes[  $f, $m , 0.5 , 32 ] -m cc[  $f, $m , 0.5 , 4 ] -t SyN[ .20, 3, 0 ] -c [ $syn ] -s 1x0.5x0vox  -f 4x2x1 -l 1 -u 1 -z 1 -o [ ${prefix} ] -v 0"


$reg -d $dim -r [ $f, $m ,1]  \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t translation[ 0.1 ] \
                         -c [$its,1.e-8,20]  \
                        -s 4x2x1vox  \
                        -f 6x4x2 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t rigid[ 0.1 ] \
                         -c [$its,1.e-8,20]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 1 , 32, regular, $percentage ] \
                         -t affine[ 0.1 ] \
                         -c [$its,1.e-8,20]  \
                        -s 4x2x1vox  \
                        -f 3x2x1 -l 1 \
                        -m mattes[  $f, $m , 0.5 , 32 ] \
                        -m cc[  $f, $m , 0.5 , 4 ] \
                         -t SyN[ .20, 3, 0 ] \
                         -c [ $syn ]  \
                        -s 1x0.5x0vox  \
                        -f 4x2x1 -l 1 -u 1 -z 1 \
                        -o [ ${prefix} ] -v 0



antsApplyTransforms -d $dim -i $m -r $f -o ${prefix}.nii -n BSpline -f 0 -v 0 -t "$prefix"0GenericAffine.mat "$prefix"1Warp.nii.gz

fslmaths $outvol -thr 0 $outvol


END=$(date +%s)
DIFF=$(( $END - $START ))
((sec=DIFF%60, DIFF/=60, min=DIFF%60, hrs=DIFF/60))
echo "${green}ANTS deformable registration took $hrs HRS $min MIN $sec SEC ${reset}"
