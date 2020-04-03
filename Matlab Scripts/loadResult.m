function pixelListCell = loadResult(filename, border, xList, yList, zList, blockID)

pixelListCell = cell(5000,1);
cellN = 0;
fidin=fopen(filename);

while ~feof(fidin)
    tline=fgetl(fidin);
    if(size(tline,2) >5)
        if(strcmp(tline(1:4),'PASS'))%%%%%
            res = tline;
            %find ':'
            id = 6;
            for i=1:size(res,2)
                if(res(i)==':')
                    id = i;
                    break;
                end
            end
            res = res(:,id+2:end);
            fidout=fopen('tmp.txt','w');
            fprintf(fidout,'%s\n\n',res);
            fclose(fidout);
            data=importdata('tmp.txt');
            data = data';
            npixels = size(data,1)/5;
            data = reshape(data,[5 npixels]);
            
            %
            pixeltmp = zeros(5,1000);
            pixelNum = 0;
            
            idx = blockID(1);
            idy = blockID(2);
            idz = blockID(3);
            
            for i=1:npixels
                x = data(3,i) - xList(idx,1) +1;
                y = data(4,i) - yList(idy,1) +1;
                z = data(5,i) - zList(idz,1) +1;
                if( (x>border)&&(x<(xList(idx,2)-border)))
                    if( (y>border)&&(y<(yList(idy,2)-border)))
                        if( (z>border)&&(z<(zList(idz,2)-border)))
                            pixelNum= pixelNum+1;
                            pixeltmp(1,pixelNum) = data(1,i);
                            pixeltmp(2,pixelNum) = data(2,i);
                            pixeltmp(3,pixelNum) = data(3,i);
                            pixeltmp(4,pixelNum) = data(4,i);
                            pixeltmp(5,pixelNum) = data(5,i);
                        end
                    end
                end
            end
            if(pixelNum >0)
                cellN = cellN + 1;
                pixelListCell{cellN,1} = pixeltmp(:,1:pixelNum);
            end
        end
    end
end
fclose(fidin);

pixelListCell = pixelListCell(1:cellN,1);
end