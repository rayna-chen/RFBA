clc
clear;

%Parameters

% load info
dirname = uigetdir('*.*','Select folder');
[file,path ] = uigetfile(  '*.mat','Select split-para file');
load(fullfile(path,file));

% Reconstruction info
zoomFactor=10;
psfSigma =3.5; %psf sigma
psfThreshold = 1000; %psf cut off threshold, relative to 65535
Output_name = [path,'\reconstruction'];

%% do Split
xList = doSplit(imgX, Block_X, Block_border);
yList = doSplit(imgY, Block_Y, Block_border);
zList = doSplit(imgZ, Block_Z, Block_border);
%% Combine

totalPixels = zeros(5,10000);
nPixel = 0;

%% Output
imgout = zeros(imgY*zoomFactor+1, imgX*zoomFactor+1, imgZ*zoomFactor+1);
psf = genpsf(psfSigma, psfThreshold);


for i=1:size(xList,1)
    for j=1:size(yList,1)
        for k=1:size(zList,1)
            
            id = [int2str(i),'_',int2str(j),'_',int2str(k)];
            resultName = [dirname, '\res_',id,'.txt'];
            blockID = zeros(3,1);
            blockID(1)= i;
            blockID(2)= j;
            blockID(3)= k;
            pixelList = loadResult(resultName, Block_border, xList, yList, zList, blockID);
            for kk=1:size(pixelList,1)
                nPixel = size(pixelList{kk,1},2);
                for ii=1:nPixel
                    px = pixelList{kk,1}(3,ii);
                    py = pixelList{kk,1}(4,ii);
                    pz = pixelList{kk,1}(5,ii);
                    
                    px1 = round(px*zoomFactor) +1;                    
                    py1 = round(py*zoomFactor) +1;                    
                    pz1 = round(pz*zoomFactor) +1;
                    
                    imgout(py1,px1,pz1) =  imgout(py1,px1,pz1)  + 1;                    
                    
                end
            end
            disp(k);
        end
    end
end
disp('combine');
%% resize imgout for gpuconv
resize_imgout=imgout(1:imgY*zoomFactor, 1:imgX*zoomFactor, 1:imgZ*zoomFactor);
Y=size(resize_imgout,1);
X=size(resize_imgout,2);
Z=size(resize_imgout,3);


%% searching neighbour
r_noise=30;
m_noise=r_noise*2+1;
zr_noise=3;
zm_noise=zr_noise*2+1;


temp_imgout=resize_imgout;
imgin_mask = double(resize_imgout>0);

zz = zr_noise+1;
xx = r_noise+1;
yy = r_noise+1;

imgVol_mask = zeros(Y,X,Z);
imgPlane_mask = zeros(Y,X,Z);
imgslice_mask = zeros(Y,X);
imgline_mask = zeros(Y,X);
tic;

for zz=1:Z
    imgslice = imgin_mask(:,:,zz);  
    for xx=1:X     
        imgline = imgslice(:,xx);
       %% Calc Line
        yy = r_noise+1;
        neighbour = imgline((yy-r_noise):(yy+r_noise));
        tempsum = sum(neighbour);
        imgline_mask(yy,xx) = tempsum;  
        for yy=(r_noise+2):(Y-r_noise)
            neighbour_plus = imgline(yy+r_noise);
            neighbour_minus = imgline(yy-r_noise-1);
            tempsum = tempsum + sum(neighbour_plus) -  sum(neighbour_minus);
            imgline_mask(yy,xx) = tempsum;
        end
    end
    
    %% Calc Slice
    xx = r_noise+1;
    neighbour = imgline_mask(:,(xx-r_noise):(xx+r_noise));
    tempsum = sum(neighbour,2);
    imgslice_mask(:,xx) = tempsum;
    for xx=(r_noise+2):(X-r_noise)
        neighbour_plus = imgline_mask(:,xx+r_noise);
        neighbour_minus = imgline_mask(:,xx-r_noise-1);
        tempsum = tempsum + neighbour_plus -  neighbour_minus;
        imgslice_mask(:,xx) = tempsum;
    end
    imgPlane_mask(:,:,zz) = imgslice_mask;
end

%% Calc Vol
zz = zr_noise+1;
neighbour = imgPlane_mask(:,:,(zz-zr_noise):(zz+zr_noise));
tempsum = sum(neighbour,3);
imgVol_mask(:,:,zz) = tempsum;

for zz=(zr_noise+2):(Z-zr_noise)
    neighbour_plus = imgPlane_mask(:,:,zz+zr_noise);
    neighbour_minus = imgPlane_mask(:,:,zz-zr_noise-1);
    tempsum = tempsum + neighbour_plus -  neighbour_minus;
    imgVol_mask(:,:,zz) = tempsum;
end

imgVol_mask2 = double(imgVol_mask >70);     % change avalue here
final_mask = ones(Y,X,Z);
final_mask((r_noise+1):(Y-r_noise),(r_noise+1):(X-r_noise),(zr_noise+1):(Z-zr_noise)) = imgVol_mask2((r_noise+1):(Y-r_noise),(r_noise+1):(X-r_noise),(zr_noise+1):(Z-zr_noise));
temp_imgout = temp_imgout.*final_mask;
clear final_mask imgin_mask imgout imgPlane_mask imgVol_mask imgVol_mask2 neighbour
toc;
disp('conv');

result = convn(temp_imgout,psf,'same');
save([Output_name,'.mat'],'temp_imgout','result');

maxres = max(max(max(result)));
result = result./maxres;
result = result.*65535;

%save tiff
for zz = 1:imgZ*zoomFactor
    ImgOut = uint16(result(:,:,zz));
    if zz==1
        imwrite(ImgOut, [Output_name,'.tif']);
    else
        imwrite(ImgOut,[Output_name,'.tif'], 'WriteMode','append');
    end
end






















