# RFBA 3D SUPER RESOLUTION

This is the software for Efficient super-resolution volumetric imaging by radial fluctuation Bayesian analysis light-sheet microscopy.

The code has a CUDA implementation for 3D images of SRRF(Super-Resolution Radial Fluctuations) gradient calculation algorithm. Part of codes are derived from NanoJ-SRRF Project. https://github.com/HenriquesLab/NanoJ-SRRF

Ref: `N. Gustafsson, S. Culley, G. Ashdown, D. M. Owen, P. M. Pereira, R. Henriques, Nat Commun 2016, 7, 12471`

The code has a multi-threaded implementation for 3D images of Bayesian analysis of fluorophore blinking and bleaching (3B) microscopy analysis. Code are base on 3B microscopy analysis software (ThreeB). https://github.com/edrosten/threeB

Ref: `aS. Cox, E. Rosten, J. Monypenny, T. Jovanovic-Talisman, D. T. Burnette, J. Lippincott-Schwartz, G. E. Jones, R. Heintzmann, Nat Methods 2011, 9(2), 195-200; bE. Rosten, G. E. Jones, S. Cox, Nat Methods 2013, 10(2), 97.`

For detailed comparison, Please check:
https://github.com/ioncannon1/threeB/commit/49aed2f7d85f47f2b3e76e9ea0fa8dddabf28c2f?diff=unified

## DIRECTORIES
RFBA Main: RFBA Main Program Source Code
Lib: Dependencies, please copy all files/folders in `/Lib` to `/RFBA Main/lib` before compliling
Matlab Scripts: MATLAB Scripts for pre-processing and post-processing
Test Data: Data for testing

## COMPILING ENVIRONMENT
Visual Studio 2015 + CUDA8.0

## BUILDING 
1. Download RFBA program(https://github.com/rayna-chen/RFBA/);
2. Copy all files/folders in `/Lib` to `/RFBA Main/lib` before compliling
3. Open srtb.sln in VS2015, Compile the solution.

## SYSTEM REQUIREMENTS
A CUDA capable NVIDIA graphics card with at least 4GiB VRAM

## TEST DATA
The test data can be download from https://github.com/rayna-chen/RFBA/tree/master/Test%20Data, which includes microtubule data and nucleus data, as well as their corresponding results.

To run the inference on your own data, make sure that:
1. The data is 3-D tiff sequence with no less than 30 volumes. 
2. The data is in 8-bit or 16-bit.



## DRIFT CORRECTION
 
Registration ImageJ-Plugin is used for drift correction(https://github.com/fiji/Descriptor_based_registration/)
1. Open multi-volume 3D-tif image in Fiji
2. Generate Hyperstacks
3. Choose Registration> Descriptor-based series registration(2d/3d+t)
* For test data, the drift has been corrected in advance.

## DATA PREPROCESSING

1. Run `GenerateBlock.m` to choose drift-corrected image and the program will prepare input data, do image split, and generate multiple batch files, which will be all saved in a new folder named by your chosen image. 
For example, if you choose nucleus.tif in the test data. the default input data is saved in `.\nucleus\img`. 
For images with size larger than 30 × 30 × 30 pixels, block computing is recommended for higher efficiency. The default block size is set as 29 × 29 × 29. The blocked markup image will be  saved in `.\nucleus\img` and the split parameters will  be saved in  `.\nucleus\Split-para.mat` for further use when combining results. 

2. Run RFBA main program
The batch file(`RFBA.bat`) includes the SRRF command at the begining and RFBA commands for every block. SRRF command need to be run at first and then you can run RFBA commands for each block parallelly. The results(`res*.txt`) will be saved in `.\img`.

3. Run `CombineBlock.m` to reconstruct image from txt files. The input includes the txt file folder(for example: `.\nucleus\img`) and `Split-para.mat` and the reconstrcution image will be saved in the parent folder of `.\img` (e.g`. .\nucleus`).

* For more detail about the parameters, please check out the annotations in MATLAB scripts.

## FAQ
For further question, please contact us feipeng@hust.edu.cn
