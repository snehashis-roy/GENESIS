===========================================================

When using this code, please cite the paper:
S. Roy, W. T. Wang, A. Carass, J. L. Prince, J. A. Butman, D. L. Pham
"PET Attenuation Correction Using Synthetic CT from Ultrashort Echo-Time MR Imaging",
Journal of Nuclear Medicine 55 (12), 2071-2077, 2014
http://jnm.snmjournals.org/content/55/12/2071.full

===========================================================


This software is intended to be used in Linux 64bit operating system

1) The main script is GENESIS.sh. It requires FSL, ROBEX, ANTs to be installed
and added to the PATH environment variable. We used ANTs 2.2.0 for all experiments.


https://fsl.fmrib.ox.ac.uk/fsl/docs/install/linux.html
http://stnava.github.io/ANTs/
http://www.nitrc.org/projects/robex


2) The script requires Matlab 2022a runtime compiler for 64bit Linux. If you do not have the
Matlab Compiler Runtime v9.12, please download it from Mathworks:

http://www.mathworks.com/products/compiler/mcr/

Unzip and install it in a suitable folder.


3) Open run_hallucination.sh, HU2umap.sh, remove_background_noise.sh,
and fix_artifacts_in_CT.sh change the MCRROOT variable to the installed
v912 directory.


4) One set of atlases is provided along with the code.
atlas_ute1       -->  1st echo of ultra-short echo time images, showing some
                      signal for bones
atlas_ute2       -->  2nd echo of ultra-short echo time images.
atlas_CT         -->  CT image of the atlas
atlas_headmask   -->  After background noise removal, binary mask of the whole head
atlas_wmmask     -->  Approximate white matter mask of the atlas, used for
                      computing median white matter intensity, which is used for
                      intensity normalization



5) Example usage: GENESIS.sh ute1.nii ute2.nii /home/user/myoutput/ 24
See GENESIS.sh for more details.

