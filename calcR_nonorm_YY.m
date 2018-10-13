% R factor calculation
% data1 and data2 should have the same size
% data1 is considered to be observed data,
% data2 is considered to be calculated data.

function R = calcR_nonorm_YY(data1,data2)

    data1 = data1(:);
    data2 = data2(:);
    
    if length(data1)~=length(data2)
        fprintf(1,'data length does not match!\n');
        R = -1;
        return
    end
    
    R = sum(abs(abs(data1)-abs(data2)))/sum(abs(data1));
    
end