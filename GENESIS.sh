#!/bin/bash

check_modules () {
    FUNC=$1    
    if ! [ -x "$(command -v "$FUNC")" ]; then
        echo "
    Error: $FUNC is not installed. Please add $FUNC to PATH.
        "
        exit 1
    fi
}
check_modules antsRegistration
check_modules antsApplyTransforms
check_modules fslmaths
check_modules remove_ext
check_modules N4BiasFieldCorrection
check_modules gzip
check_modules runROBEX.sh
check_modules ROBEX
check_modules N4BiasFieldCorrection



if [ $# -lt "3" ]; then
    echo "----------------------------------------------------------------
    Usage:
    ./GENESIS.sh UTE1 UTE2 OUTPUTDIR NUMCPU
    UTE1            Subject UTE 1st echo (visible bone signal) NIFTI .nii file
    UTE2            Subject UTE 2nd echo (no bone signal) NIFTI .nii file
    OUTPUTDIR       Output directory where final results and all temporary
                    files will be written
    NUMCPU          Number of parallel processing cores to use. Maximum 12

    All files must be NIFTI .nii
    This script needs FSL, ROBEX, ANTS to be installed and added to PATH.
    
    ----------------------------------------------------------------"
    exit 1
fi



export AFNI_NIFTI_TYPE_WARN=NO
ROOT="$0"
ROOT=`dirname $ROOT`
ROOT=`readlink -f $ROOT`
UTE1=$1
UTE2=$2
OUTDIR=$3
NUMCPU=$4


if [ x"$NUMCPU" == "x" ];then
    NUMCPU=12
fi
echo "Using $NUMCPU parallel processing cores"
export ITK_GLOBAL_DEFAULT_NUMBER_OF_THREADS=$NUMCPU
export FSLOUTPUTTYPE=NIFTI
START=$(date +%s)

UTE1=`readlink -f $UTE1`
UTE2=`readlink -f $UTE2`
OUTDIR=`readlink -f $OUTDIR`
red=`tput setaf 5`
green=`tput setaf 2`
reset=`tput sgr0`


ATLASDIR=$ROOT/atlas/
ATLASWMMASK=$ATLASDIR/atlas_wmmask.nii
ATLASHMASK=$ATLASDIR/atlas_headmask.nii


ID=`basename $UTE1`
ID=`remove_ext ${ID}`


mkdir -p ${OUTDIR}
if [ ! -d $OUTDIR/${ID} ];then
    mkdir $OUTDIR/${ID}
else
    echo "${red}WARNINH: Working directory $OUTDIR/${ID} already exists, its content may be overwritten."
    echo "If necessary, please clean the existing directoy, and restart. ${reset}"
    sleep 5
fi

echo "Working  directory: ${green}$OUTDIR/${ID}/ ${reset}"
cd $OUTDIR/${ID}/


#ORI=`3dinfo -orient $UTE1`
#3dresample -orient RAI -inset $UTE1 -prefix $OUTDIR/${ID}/ute1.nii
#3dresample -orient RAI -inset $UTE2 -prefix $OUTDIR/${ID}/ute2.nii
fslmaths $UTE1  $OUTDIR/${ID}/ute1.nii -odt float
fslmaths $UTE2 $OUTDIR/${ID}/ute2.nii -odt float
echo N4BiasFieldCorrection -s 3 -d 3 -c [ 50x50x50x50,0.00001] -i ute1.nii -o ute1.nii -v 0 -r 1
N4BiasFieldCorrection -s 3 -d 3 -c [ 50x50x50x50,0.00001] -i ute1.nii -o ute1.nii -v 0  -r 1
echo N4BiasFieldCorrection -s 3 -d 3 -c [ 50x50x50x50,0.00001] -i ute2.nii -o ute2.nii -v 0 -r 1
N4BiasFieldCorrection -s 3 -d 3 -c [ 50x50x50x50,0.00001] -i ute2.nii -o ute2.nii -v 0 -r 1

echo "Register ute2 to ute1 via rigid registration --"
$ROOT/antsaffine.sh ute1.nii ute2.nii ute2_reg.nii yes
mv -vf ute2_reg.nii ute2.nii


#echo "Register atlas to ute2 via approximate ANTS to remove background noise and estimate white matter intensity --"
#$ROOT/AntsExample.sh ute2.nii $ATLASDIR/atlas_ute2.nii fastfortesting atlas_reg.nii
#antsApplyTransforms -d 3 -i $ATLASWMMASK -r ute2.nii -o atlas_wmmask_reg.nii -n NearestNeighbor -t atlas_reg0GenericAffine.mat atlas_reg1Warp.nii.gz  -f 0 -v 0
#antsApplyTransforms -d 3 -i $ATLASHMASK -r ute2.nii -o atlas_headmask_reg.nii -n NearestNeighbor -t atlas_reg0GenericAffine.mat atlas_reg1Warp.nii.gz  -f 0 -v 0
#rm -f *Warp.nii.gz

echo "Skullstrip ute2 to find robust normalization factor --"
echo runROBEX.sh ute2.nii ute2_strip.nii   # skullstrip ute2 and use brainmask for fixing some artifacts after synthesis
runROBEX.sh ute2.nii ute2_strip.nii 
echo fslmaths ute2_strip.nii -bin brainmask.nii
fslmaths ute2_strip.nii -bin brainmask.nii
gzip -vf ute2_strip.nii

#echo "Dilating headmask for robust removal of background noise --"
#for u in 1 2 3 4
#do
#    fslmaths atlas_headmask_reg.nii -kernel boxv 3 -dilM atlas_headmask_reg.nii -odt char
#done
#echo fslmaths ute1.nii -mul atlas_headmask_reg.nii ute1.nii
#fslmaths ute1.nii -mul atlas_headmask_reg.nii ute1.nii


# Removal of background noise is crucial because synthesizing the noise voxels could take long time
echo "${green}$ROOT/remove_background_noise.sh ute1.nii 30 ute1_denoised.nii ${reset}"
$ROOT/remove_background_noise.sh ute1.nii 30 ute1_denoised.nii


# create noise removed mask and images
fslmaths ute1_denoised.nii -bin ute1_mask.nii -odt char
fslmaths ute2.nii -mul ute1_mask.nii ute2_denoised.nii


##########
# Find mean of the images based on wmmask of atlas
#fslmaths ute2_denoised.nii -mul atlas_wmmask_reg.nii test.nii
#p2=`fslstats test.nii -P 50`
#rm -f test.nii
#echo "Computing WM median intensities : "
#fslmaths ute1_denoised.nii -mul atlas_wmmask_reg.nii test.nii
#p1=`fslstats test.nii -P 50`
#rm -f test.nii
#p1=`echo $p1 | tr -d ' '`
#p2=`echo $p2 | tr -d ' '`
##########

p1=`$ROOT/find_WM_peak.sh ute1_denoised.nii t1 brainmask.nii`
p2=`$ROOT/find_WM_peak.sh ute2_denoised.nii t1 brainmask.nii`

echo "WM peak intensity for UTE1 is found at ${p1},"
echo "WM peak intesnity for UTE2 is found at ${p2},"


# synthesis
echo "${green}run_hallucination.sh $ATLASDIR/atlas_ute1.nii,$ATLASDIR/atlas_ute2.nii,$ATLASDIR/atlas_CT.nii ute1_denoised.nii,ute2_denoised.nii no 6E4 1E5 12 3x3x1 $NUMCPU none T1,T1 CT 115 0 88.50,86.01 $p1,$p2 "${ID}"_synth.nii ${reset}"
$ROOT/run_hallucination.sh $ATLASDIR/atlas_ute1.nii,$ATLASDIR/atlas_ute2.nii,$ATLASDIR/atlas_CT.nii ute1_denoised.nii,ute2_denoised.nii no 6E4 1E5 12 3x3x1 $NUMCPU none T1,T1 CT 115 0 88.50,86.01 $p1,$p2 "${ID}"_synth.nii

echo "Smoothing with a 0.5 mm Gaussian kernel"
fslmaths "${ID}"_synth.nii -s 0.50 "${ID}"_synth_blur.nii



# Fix small artifacts in CT
echo fix_artifacts_in_CT.sh "${ID}"_synth_blur.nii brainmask.nii "${ID}"_synth_blur.nii
$ROOT/fix_artifacts_in_CT.sh "${ID}"_synth_blur.nii brainmask.nii "${ID}"_synth_blur.nii


# create mu-map from CT
echo $ROOT/HU2umap.sh "${ID}"_synth_blur.nii "${ID}"_synth_mumap.nii  120
$ROOT/HU2umap.sh "${ID}"_synth_blur.nii "${ID}"_synth_mumap.nii  120



fslmaths "${ID}"_synth_mumap.nii -mul ute1_mask.nii "${ID}"_synth_mumap.nii
fslmaths "${ID}"_synth_mumap.nii -thr 0 "${ID}"_synth_mumap.nii


gzip -vf *.nii
echo "Working  directory: ${green}$OUTDIR/${ID}/ ${reset}"
END=$(date +%s)
DIFF=$(( $END - $START ))

((sec=DIFF%60, DIFF/=60, min=DIFF%60, hrs=DIFF/60))
if [ $hrs != "0" ];then
    echo "${green}GENESIS CT Synthesis took $hrs hours $min minutes and $sec seconds ${reset}"
else
    echo "${green}GENESIS CT Synthesis took $min minutes and $sec seconds ${reset}"
fi



