



<!-- PROJECT LOGO -->
<br />
<div align="center">
  

  <h3 align="center">GENESIS - Generative Sub-Image Synthesis</h3>

  
</div>



<!-- TABLE OF CONTENTS -->
<details>
  <summary>Table of Contents</summary>
  <ol>
    <li>
      <a href="#about-the-project">About The Project</a>      
    </li>
    <li>
      <a href="#getting-started">Getting Started</a>
      <ul>
        <li><a href="#prerequisites">Prerequisites</a></li>
        <li><a href="#installation">Installation</a></li>
      </ul>
    </li>
    <li><a href="#usage">Usage</a></li>
    <li><a href="#reference">Reference</a></li>       
    <li><a href="#contact">Contact</a></li>
    
  </ol>
</details>



<!-- ABOUT THE PROJECT -->
## About The Project
<p align="center">
  <img src="https://github.com/snehashis-roy/GENESIS/blob/master/F1.jpg"  height="800"/>  
</p>


GENESIS (**GEN**erative **S**ub-**I**mage **S**ynthesis) is a patch (or sub-image) based
algorithm to learn intensity mapping between dual ultrashort echo-time (UTE) images and CT. 

When using this code, please cite the paper:
```
S. Roy, W. T. Wang, A. Carass, J. L. Prince, J. A. Butman, D. L. Pham
"PET Attenuation Correction Using Synthetic CT from Ultrashort Echo-Time MR Imaging",
Journal of Nuclear Medicine 55 (12), 2071-2077, 2014
http://jnm.snmjournals.org/content/55/12/2071.full
```

The algorithm uses an atlas and a subject image set. The atlas image set comprises of a pair of UTE and 
the corresponding registered CT images. A subject image set contains a pair of UTE images. 
A transformation is learned between a subject patch cloud and an atlas patch cloud to 
generate a *synthetic*  CT patch of the subject.

While a *synthetic* CT image is not intended for clinical diagnosis, it can be used for 
various image analysis tasks, such as registration or attenuation correction. See [1,2,3]
for applications.




<!-- GETTING STARTED -->
## Getting Started
Although all Matlab source files are provided, GENESIS is intended to be used as
a command line tool via a Bash script ```GENESIS.sh``` that takes two UTE subject images
and runs the complete pipeline, assuming the atlas is located in the same folder as 
the script. The pipeline comprises of compiled codes using Matlab R2022a compiler on 64bit Linux.

There are 3 main steps of the pipeline,
1. N4 bias correction and affine registration between UTE images. A deformable registration
   script (```AntsExample.sh```) is also provided if affine is not sufficient.
3. Computing robust normalization factors by skull-stripping one of the input channels and finding
   the WM peak intensity (i.e., mode). See ```find_WM_peak.sh``` script.
4. Use the two atlas channels along with their normalization factors as well as
   the subject channels to create a patch cloud and compute synthetic CT intensities
   of the subject. See ```run_hallucination.sh``` and ```image_hallucination.m```
   for details.

### Prerequisites

Please install FSL, ROBEX, MCR, and ANTs and add the binaries to the PATH variable

* FSL : 
  ```https://fsl.fmrib.ox.ac.uk/fsl/docs/install/linux.html```
* ANTs : ```https://github.com/ANTsX/ANTs/releases```
  
We have extensively tested the pipeline with ANTs version 2.2.0. Note that with
newer versions of ANTs, it is possible that some command line arguments
we used may be deprecated.

* ROBEX :  ```http://www.nitrc.org/projects/robex```
  
* Matlab Compiler Runtime (MCR) : Download the 64-bit Linux MCR installer for MATLAB 2022a (v912)
```
https://www.mathworks.com/products/compiler/matlab-runtime.html
```

FSL's fslmaths command is used for various image processing tasks like
image multiplication etc. ANTs registration is used to register multiple
modalities in an affine manner. ROBEX is used to skullstrip the input 
images to find WM peak intensity (i.e., mode) as a robust normalization factor.



### Installation
1. Install the MCR (v912) to somewhere suitable.

2. Add the MCR installation path, i.e. the v912 directory, to the included compiled shell scripts' MCRROOT variable. 
In each of the following 4 shell scripts, remove_background_noise.sh, fix_artifacts_in_CT.sh, HU2umap.sh, run_hallucination.sh,
replace the line containing ```MCRROOT=/usr/local/matlab-compiler/v912```
to the path where the MCR is installed, e.g.,
```
MCRROOT=/home/user/MCR/v912
```
if the MCR is installed in ```/home/user/MCR```.

3. Install ANTs binaries and add the binary path to shell \$PATH
```
export PATH=/home/user/ANTs-2.2.0/install/bin:${PATH}
```

4. Download the latest ROBEX binaries (ROBEXv12.linux64.tar.gz), and unzip to
   suitable folder, then add that to \$PATH.
5. Install FSL either using ```curl``` or via Python
   ```python fslinstaller.py -d /home/user/FSL/```
   Then add the corresponding bin folder to \$PATH
   ```export PATH=${PATH}:/home/user/FSL/bin```





<!-- USAGE EXAMPLES -->
## Usage

The main script is ```GENESIS.sh```. If all the dependencies are satisfied, ```./GENESIS.sh```
should show the usage:
```
 Usage:
    ./GENESIS.sh UTE1 UTE2 OUTPUTDIR NUMCPU
    UTE1            Subject UTE 1st echo (visible bone signal) NIFTI .nii file
    UTE2            Subject UTE 2nd echo (no bone signal) NIFTI .nii file
    OUTPUTDIR       Output directory where final results and all temporary
                    files will be written
    NUMCPU          Number of parallel processing cores to use. Maximum 12

    All files must be NIFTI .nii
    This script needs FSL, ROBEX, ANTS to be installed and added to PATH.
```



### Notes:
1. The ```GENESIS.sh``` script is hardcoded to accept 2 modalities. However,
the main script ```run_hallucination.sh``` can take any number of modalities.
2. If using more than 2 modalities to synthesize CT, please run the
   ```run_hallucination.sh``` script using appropriate arguments. E.g.,
   assuming T1,T2, and FLAIR images are used to synthesize a CT,
   
```
./run_hallucination.sh atlas_t1.nii,atlas_t2.nii,atlas_fl.nii,atlas_CT.nii  subjt1.nii,subjt2.nii,subjfl.nii \
      no 6E4 1E5 12 3x3x1 $NUMCPU none \
      T1,T2,FL CT 115 0  $Ap_1,$Ap_2,$Ap_3 $Sp_1,$Sp_2,$Sp_3  subjct_synth.nii
```
where Ap<sub>x</sub> and Sp<sub>x</sub>, ```x=1,2,3``` denote the atlas and subject WM peak intensities
that can be found using ```find_WM_peak.sh``` script.


<!-- Reference -->
## Reference
1. S. Roy, A. Carass, A. Jog, et al., "MR to CT registration of brains using image synthesis", Proc SPIE. March 21, 2014:9034. ([Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC4104818/))
2. S. Roy, A. Jog, A. Carass, J. L. Prince, "Atlas based intensity transformation of brain MR images", Multimodal Brain Image Analysis. 2013;8159:51–62. ([Link](https://link.springer.com/chapter/10.1007/978-3-319-02126-3_6))
3. S. Roy, Y.-Y. Chou, A. Jog, J. A. Butman, D. L. Pham, "Patch Based Synthesis of Whole Head MR Images: Application To EPI Distortion Correction",
   Intl. Workshop on Simulation and Synthesis in Medical Imaging. 2016;9968:146–156. ([Link](https://pmc.ncbi.nlm.nih.gov/articles/PMC5375117/))


<!-- CONTACT -->
## Contact

Snehashis Roy - snehashis[dot]roy@gmail[dot]com

Project Link: [https://github.com/snehashis-roy/GENESIS](https://github.com/snehashis-roy/GENESIS)


<p align="right">(<a href="#readme-top">back to top</a>)</p>



