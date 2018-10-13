% origin of model should be (0,0,0)
% model is 3xN array, with real unit (A, nm, etc.)
% Res is pixel size, same unit with model

function [cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D(model, atomtype, currAtomInd, projAngles, Heights, Bfactors, totVolSize,Res, CropHalfWidth, cylinderRad, cylinderAtomBufferRad)

model = model ./ Res;
%FHeights = Heights;
%FWidths = Bfactors / pi^2 / Res^2;

%[tX,tY] = ndgrid(-cylinderRad:cylinderRad,-cylinderRad:cylinderRad);
%circInd = find( tX.^2+tY.^2 < (cylinderRad+0.5)^2);


if length(totVolSize) == 1
  totVolSize = [totVolSize totVolSize totVolSize];
end
  
  
cenPos = round( (totVolSize +1)/2);
cylinderProjAr = zeros(cylinderRad*2+1,cylinderRad*2+1,size(projAngles,1));
cylinerIndsAr = zeros(2,cylinderRad*2+1,size(projAngles,1));

for i=1:size(projAngles,1)
  currAngle = projAngles(i,:);
  
  RotMat1 = MatrixQuaternionRot([0 0 1],currAngle(1));
  RotMat2 = MatrixQuaternionRot([0 1 0],currAngle(2));
  RotMat3 = MatrixQuaternionRot([1 0 0],currAngle(3));
  
  
  
  ROTmodel = (RotMat1*RotMat2*RotMat3)'*model;
  
  currAtomPos= ROTmodel(:,currAtomInd);
  
  %currAtomXY = currAtomPos(1:2);
  currAtomXY_round = round(currAtomPos(1:2));
  
  cropIndX = (-cylinderRad:cylinderRad) + currAtomXY_round(1) + cenPos(1);
  cropIndY = (-cylinderRad:cylinderRad) + currAtomXY_round(2) + cenPos(2); 
  
  cylinderAtomInd = find( (ROTmodel(1,:)-currAtomXY_round(1)).^2 + (ROTmodel(2,:)-currAtomXY_round(2)).^2 < (cylinderAtomBufferRad+cylinderRad+0.5)^2);
  %cylinderAtomBufferRad
  %circInd
  
  inputModel =ROTmodel(:,cylinderAtomInd) - repmat([currAtomXY_round(1);currAtomXY_round(2);0],[1 length(cylinderAtomInd)]);
  inputAtomType = atomtype(cylinderAtomInd);
  
  Vol = My_create_volProj_from_model_exact_HB_fixedfa_2D(inputModel*Res, inputAtomType, Heights, Bfactors, [cylinderRad*2+1 cylinderRad*2+1 totVolSize(3)] + cylinderRad*2, Res, CropHalfWidth);
  
  Vol = My_stripzero(Vol,[cylinderRad*2+1 cylinderRad*2+1]);
  cylinderProjAr(:,:,i) = Vol;
  cylinerIndsAr(:,:,i) = [cropIndX;cropIndY];
end
  
