clear;
clc;

%% Parameters needed here 

% Input image infromation
Z_slice=15;
Vol=30;

% Block size
Block_X = 23;
Block_Y = 23;
Block_Z = 15;
Block_border = 3;


% Algroithm input
maxIter = 50;  %Maximum Iteration
zfactor =0.5764; %Axial pixel number occupied by one PSF/ lateral pixel number occupied by one PSF
pixelBasedInitialPoints = 0;    % using Pixel Num based Initial Points, or using Intensity based Initial Points
initPointsFactor = 38.9*4;      %0--initial num of points = pixel num of block / initPointsFactor
                                %1--initial num of points =pixel num of block / initPointsFactor * initPointsFactor_int


%% Cut hyperstack into multiple tiff image

[filename, filepath ] = uigetfile(  '*.tif','Select  images');
info     = imfinfo(fullfile(filepath, filename));
N_slice = size(info,1);
height   = info.Height;
width    = info.Width;
bitdepth = info.BitDepth;
assert (bitdepth == 8 || bitdepth == 16, 'Only 8bit and 16bit are supported')
assert (Z_slice * Vol == N_slice, 'Slice number dose not match')
save_dir = [filename(1:(length(filename)-4)),'\img'];
if ~exist(save_dir, 'dir')
    mkdir(save_dir)
end

for i=1:N_slice
    
   img(:,:,i)=imread(fullfile(filepath,filename),i);
end

for i=1:Vol
  cut_img(:,:,:)=img(:,:,(i-1)*Z_slice+1:(i-1)*Z_slice+Z_slice);      
   
  imwrite(cut_img(:,:,1),[save_dir,'\cell',num2str(i),'.tiff']);    
 for k=2:Z_slice
    imwrite(cut_img(:,:,k),[save_dir,'\cell',num2str(i),'.tiff'],'WriteMode','append');   
 end
end 

%% Generate Block mask and bat file

[imgY,imgX,imgZ]=size(cut_img);
dirname=save_dir;
imgfiles = [dirname,'\cell*'];
avgimage =  [dirname,'\cell1.tiff'];  %average of image, or just 1 image


% Load Average Image
avgImg = zeros(imgY,imgX, imgZ);
for i=1:imgZ
    avgImg(:,:,i) = imread(avgimage, i);
end

% Split
xList = doSplit(imgX, Block_X, Block_border);
yList = doSplit(imgY, Block_Y, Block_border);
zList = doSplit(imgZ, Block_Z, Block_border);


%% Save split information
name= [filename(1:(length(filename)-4)),'\Split-para.mat'];
save( name, 'Block_X', 'Block_Y' ,'Block_Z', 'Block_border', 'imgX' ,'imgY' ,'imgZ');

%% Gen Markup and cmd
fp = fopen('RFBA.bat','w');
resultName = [dirname, '\resdummy.txt'];
markupName = [dirname, '\markupdummy.tif'];
for zz = 1:imgZ
    markupImgOut = uint8(ones(imgY, imgX, 3));
    markupImgOut = markupImgOut.*255;
    if zz==1
        imwrite(markupImgOut, markupName);
    else
        imwrite(markupImgOut, markupName, 'WriteMode','append');
    end
end
cmdline = ['@.\RFBA --srrf_only 1 --save_spots ', resultName, '  --log_ratios ', markupName, ' ', imgfiles];
fprintf(fp, '%s\n', cmdline);

%% calculate initPointsFactor_int

for i=1:size(xList,1)
    for j=1:size(yList, 1)
        for k=1:size(zList,1)
            
            
            x1 = xList(i,1);
            x2 = x1 + xList(i,2) -1;
            y1 = yList(j,1);
            y2 = y1 + yList(j,2) -1;
            z1 = zList(k,1);
            z2 = z1 + zList(k,2) -1;
            sum_intensity(i,j,k)= sum(sum(sum(avgImg(y1:y2,x1:x2,z1:z2))));
             
        end
    end
end

  av_intensity=sum_intensity/(xList(i,2)* yList(j,2)*zList(k,2));
  normalized_inensity=av_intensity/(mean(mean(mean(av_intensity))));
  initPointsFactor_int=sqrt(normalized_inensity);
  
%% generate mask and bat  
for i=1:size(xList,1)
    for j=1:size(yList, 1)
        for k=1:size(zList,1)
            id = [int2str(i),'_',int2str(j),'_',int2str(k)];
            resultName = [dirname, '\res_',id,'.txt'];
            markupName = [dirname, '\markup_',id,'.tif'];
            x1 = xList(i,1);
            x2 = x1 + xList(i,2) -1;
            y1 = yList(j,1);
            y2 = y1 + yList(j,2) -1;
            z1 = zList(k,1);
            z2 = z1 + zList(k,2) -1;
            markupImg = uint8(zeros(imgY, imgX, imgZ));
            markupImg(y1:y2,x1:x2, z1:z2) = 255;
            
            for zz = 1:imgZ
                markupImgOut = uint8(zeros(imgY, imgX, 3));
                markupImgOut(:,:,1)= markupImg(:,:,zz);
                markupImgOut(:,:,2)= markupImg(:,:,zz);
                markupImgOut(:,:,3)= markupImg(:,:,zz);
                if zz==1
                    imwrite(markupImgOut, markupName);
                else
                    imwrite(markupImgOut, markupName, 'WriteMode','append');
                end
            end
            

            
            
           
            iniPoints = ceil(xList(i,2)* yList(j,2)*zList(k,2)/initPointsFactor);
                
            if(pixelBasedInitialPoints==0)  
            iniPoints = ceil(iniPoints*initPointsFactor_int(i,j,k));
            end

            cmdline = ['@.\RFBA --z_factor ', num2str(zfactor), ' --split_mode 1 --add_guide 1 --max_iter ' ,int2str(maxIter),' --start_points ',int2str(iniPoints), ' '];
            cmdline = [cmdline, '--save_spots ', resultName, '  --log_ratios ', markupName, ' ', imgfiles];
            
            fprintf(fp, '%s\n', cmdline);
            
        end
    end
end

fprintf(fp, '@pause;\n;\n');

fclose(fp);