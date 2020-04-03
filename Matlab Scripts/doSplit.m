function posList  = doSplit(xall, maxX, border)

if(maxX>=xall)
    posList = [1 xall]; 
else
    %do split
    posList = [1 maxX]; 
    nowID = 1;
    while (posList(nowID) + maxX)<=xall
        nowID = nowID + 1;
        posList(nowID, 1) =  posList(nowID-1, 1) + posList(nowID-1, 2) - 2*border;   
        posList(nowID, 2) = maxX;
    end
    posList(nowID, 2) = xall -  posList(nowID, 1) + 1; 
    
end



end