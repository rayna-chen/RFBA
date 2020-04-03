# RFBA 3D super resolution

Efficient super-resolution volumetric imaging by radial fluctuation Bayesian analysis light-sheet microscopy

Codes are modefied based on 3B algorithm and SRRF algorithm. For detailed comparison, Please check out https://github.com/ioncannon1/threeB/commit/49aed2f7d85f47f2b3e76e9ea0fa8dddabf28c2f?diff=unified


## COMPILING ENVIRONMENT

Visual Studio 2015+CUDA8.0

## DATA REQUIREMENT
The test data can be download from XXXX, which includes simulated microtubule data and experimental nucleus data.
To run the inference on your own data, make sure that:
1. The data is 3-D tiff sequence with no less than 30 volumes. 
2. The data is in 8-bit or 16-bit.

## MANUAL COMPILATION 

1. Download RFBA program(https://github.com/rayna-chen/RFBA/);
2. Copy all files in Lib folder to the lib folder beneath RFBA directory.
3. Compile RFBA program.


## DRIFT CORRECTION
 
Registration ImageJ-Plugin is used for drift correction(https://github.com/fiji/Descriptor_based_registration/)
1. Open multi-volume 3D-tif image in Fiji
2. Generate Hyperstacks
3. Choose Registration> Descriptor-based series registration(2d/3d+t)
* For test data, the drift has been corrected in advance.

## DATA PREPROCESSING

1. Run GenerateBlock.m to choose drift-corrected image and the programe will prepare input data, do image split, and generate multiple batch files, which will be all sved in new folder named by your chosen image. 
For example, if you choose nucleus.tif in the test data. the default input data is saved in .\nucleus\img. 
For images with size larger than 30 × 30 × 30 pixels, block computing is recommended for higher efficiency. The default block size is set as 29 × 29 × 29. The blocked markup image will be  saved in .\nucleus\img and the split parameters will  be saved in  .\nucleus\Split-para.mat for further use when combining results. 

2. Run RFBA main program
The batch file(RFBA.bat) includes the SRRF commond at the begining and RFBA commonds for every block. SRRF commond need to be run at first and then you can run RFBA commonds for each block parallelly. The results(res*.txt) will be saved in .\img.

3. Run Combine.m to reconstruct image from txt files. The input includes the txt file folder(for example: .\nucleus\img) and Split-para.mat and the reconstrcution image will be saved in the parent folder of .\img (e.g. .\nucleus).

* For more details about the parameter, please check out the annotations in matlab scripts.

## FAQ
For further question, please contact us feipeng@hust.edu.cn

