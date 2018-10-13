% origin of model should be (0,0,0)
% model is 3xN array, with real unit (A, nm, etc.)
% Res is pixel size, same unit with model

function Vol = My_create_volProj_from_model_exact_HB_fixedfa_2D(model, atomtype, Heights, Bfactors, volsize, Res, CropHalfWidth, fixedfa)

model = model ./ Res;
FHeights = Heights;
FWidths = Bfactors / pi^2 / Res^2;

if length(volsize) == 3
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = (1:volsize(2)) - round((volsize(2)+1)/2);
    z = (1:volsize(3)) - round((volsize(3)+1)/2);
elseif length(volsize) == 1
    x = (1:volsize(1)) - round((volsize(1)+1)/2);
    y = x;
    z = x;
else
    error('volsize should be either length 3 or length 1!')
end
    
sizeX = [length(x) length(y)];

inInd = find(model(1,:) >= min(x) & model(1,:) <= max(x) & ...
             model(2,:) >= min(y) & model(2,:) <= max(y) & ...
             model(3,:) >= min(z) & model(3,:) <= max(z));
         
calcModel = model(:,inInd);
calcAtomtype = atomtype(:,inInd);

finalProj_padded = zeros( [sizeX + (CropHalfWidth+1)*2, length(Heights)]);

cenPos = round((size(finalProj_padded)+1)/2);
cropIndRef = -CropHalfWidth:CropHalfWidth;
[cropX,cropY] = ndgrid(cropIndRef,cropIndRef);

currPos_all = calcModel(1:2,:) + repmat(cenPos(1:2)',[1, size(calcModel,2)]);
currRndPos_all = round(currPos_all);
diffPos_all = currPos_all - currRndPos_all;
diffPosZ_all =  calcModel(3, :) - round(calcModel(3, :));


% try
%     pool = gcp();
%     num_threads = pool.NumWorkers;
% catch
%     num_threads = 4;
%     pool = parpool(num_threads);
% end

finalProj_padded = zeros( [sizeX + (CropHalfWidth+1)*2, length(Heights)]);
% finalProjs_padded = cell(1, num_threads);
% for t = 1:num_threads
%     finalProjs_padded{t} = finalProj_padded;
% end
% finalProjs_padded = cell(1, size(calcModel,2));
% for t = 1:size(calcModel,2)
%     finalProjs_padded{t} = finalProj_padded;
% end
% for i=1:size(calcModel,2)
%     which_thread_index = 1 + mod(i, num_threads);
%     currPos = currPos_all(:, i);
%     currRndPos = currRndPos_all(:, i);
% %     currPos = calcModel(1:2,i) + cenPos(1:2)';
%     currRndPos = round(currPos);
%     
%     cropInd1 = cropIndRef + currRndPos(1);
%     cropInd2 = cropIndRef + currRndPos(2);
%     
% %     CropProj = finalProj_padded(cropInd1,cropInd2,calcAtomtype(i));
%     diffPos = diffPos_all(:, i);
%     diffPosZ = diffPosZ_all(i);
% %     diffPos = currPos-currRndPos;
% %     diffPosZ = calcModel(3,i) - round(calcModel(3,i));
%     
%     gaussCalc = FHeights(calcAtomtype(i))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 )/FWidths(calcAtomtype(i)) );
%     
%     gaussZcalc = exp(-1*(cropIndRef - diffPosZ).^2 / FWidths(calcAtomtype(i)) );
% %     finalProj_padded(cropInd1,cropInd2,calcAtomtype(i)) = CropProj + gaussCalc*sum(gaussZcalc);
% %     finalProjs_padded{which_thread_index}(cropInd1,cropInd2,calcAtomtype(i)) = finalProj_padded{which_thread_index}(cropInd1,cropInd2,calcAtomtype(i)) + gaussCalc*sum(gaussZcalc);
% % finalProjs_padded{which_thread_index} = finalProjs_padded{which_thread_index} + 1;
% %     finalProjs_padded{i} = finalProjs_padded{i} + 1;
%     finalProjs_padded{i}(cropInd1,cropInd2,calcAtomtype(i)) = finalProjs_padded{i}(cropInd1,cropInd2,calcAtomtype(i)) + gaussCalc*sum(gaussZcalc);
% 
% end
% 
% for t = 1:size(calcModel,2)
% %     finalProj_padded = sum(finalProjs_padded,4);
%     finalProj_padded = finalProj_padded + finalProjs_padded{t};
% end

for i=1:size(calcModel,2)
    
%     currPos = currPos_all(:, i);
    currRndPos = currRndPos_all(:, i);
%     currPos = calcModel(1:2,i) + cenPos(1:2)';
%     currRndPos = round(currPos);
    
    cropInd1 = cropIndRef + currRndPos(1);
    cropInd2 = cropIndRef + currRndPos(2);
    
%     CropProj = finalProj_padded(cropInd1,cropInd2,calcAtomtype(i));
    diffPos = diffPos_all(:, i);
    diffPosZ = diffPosZ_all(i);
%     diffPos = currPos-currRndPos;
%     diffPosZ = calcModel(3,i) - round(calcModel(3,i));
    
    gaussCalc = FHeights(calcAtomtype(i))*exp( -1*( (cropX-diffPos(1)).^2 + (cropY-diffPos(2)).^2 )/FWidths(calcAtomtype(i)) );
    
    gaussZcalc = exp(-1*(cropIndRef - diffPosZ).^2 / FWidths(calcAtomtype(i)) );
%     finalProj_padded(cropInd1,cropInd2,calcAtomtype(i)) = CropProj + gaussCalc*sum(gaussZcalc);
    finalProj_padded(cropInd1,cropInd2,calcAtomtype(i)) = finalProj_padded(cropInd1,cropInd2,calcAtomtype(i)) + gaussCalc*sum(gaussZcalc);

    
end

finalProj_summed = zeros( sizeX);
if nargin < 8
    kx = 1:size(finalProj_summed,1);
    ky = 1:size(finalProj_summed,2);

    MultF_X = 1/(length(kx)*Res);
    MultF_Y = 1/(length(ky)*Res);

    CentPos = round((size(finalProj_summed)+1)/2);
    [KX, KY] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y);
    q2 = KX.^2 + KY.^2;
    clear KX KY KZ

    fa78 = fatom_vector( sqrt(q2),78 );
    fa26 = fatom_vector( sqrt(q2),26 );
    % fa26(1)
    % fa78(1)
    fixedfa = 0.5*(fa78+fa26);
end
for j=1:length(Heights)
    %fa = reshape(fatom_vector( sqrt(q2),AtomNumbers(j)),sizeX);
    CVol = My_stripzero_ajp(finalProj_padded(:,:,j),sizeX);
    FVol = My_FFTN(CVol);
    FVol = FVol .*  reshape(fixedfa,sizeX) ;
    finalProj_summed =finalProj_summed+FVol;
end

Vol = real(My_IFFTN(finalProj_summed));
end