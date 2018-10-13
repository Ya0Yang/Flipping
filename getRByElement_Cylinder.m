function R_element = getRByElement_Cylinder(model, curr_atoms_INC, curr_IncInd_i, element_flag, angles, htAr, bfAr, volSize, Res, CropHalfWidth, cylinderRad, cylinderAtomBufferRad, fixedfa_small, projections, calcR_halfRad, useInd, SF)
curr_atoms_INC(curr_IncInd_i) = element_flag;

%         [cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D(model, curr_atoms_INC, curr_IncInd(i), angles, [htAr],[bfAr ],volSize, Res, CropHalfWidth, round(cylinderRad/0.338), round(cylinderAtomBufferRad/0.338), make_fixedfa(volSize, Res));
[cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D(model, curr_atoms_INC, curr_IncInd_i, angles, [htAr],[bfAr ],volSize, Res, CropHalfWidth, round(cylinderRad/0.338), round(cylinderAtomBufferRad/0.338), fixedfa_small);

subProjs = zeros(size(cylinerIndsAr,2),size(cylinerIndsAr,2),size(cylinerIndsAr,3));
for j=1:size(angles,1)
  subProjs(:,:,j) = projections(cylinerIndsAr(1,:,j),cylinerIndsAr(2,:,j),j);
end

subProjs_crop = My_stripzero_ajp(subProjs ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
cylinderProjAr_crop = My_stripzero_ajp(cylinderProjAr ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);


R_element = calcR_nonorm_YY(subProjs_crop(useInd),cylinderProjAr_crop(useInd)*SF);