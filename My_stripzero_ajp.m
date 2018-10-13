function Vol = My_stripzero_ajp(PadVol,orisize)
    if length(size(PadVol))~=length(orisize)
        fprintf(1,'input volume dimension and and paddedsize length does not match!\n');
    elseif sum(size(PadVol)<orisize) > 0
        fprintf(1,'paddedsize should be equal to or smaller than original volume in all dimensions!\n');
    else  
nc = round((orisize+1)/2);
nc_vol = round((size(PadVol)+1)/2);
vecX = (1:orisize(1)) - nc(1) + nc_vol(1);
if numel(nc) > 1
    vecY = (1:orisize(2)) - nc(2) + nc_vol(2);
end
if numel(nc) > 2
    vecZ = (1:orisize(3)) - nc(3) + nc_vol(3);
end   
if numel(nc) == 2
    Vol = PadVol(vecX, vecY);
else
    Vol = PadVol(vecX, vecY, vecZ);
end
end