function R_elements = par_getRByElement_Cylinder(model, curr_atoms_INC, curr_IncInd_i, element_flags, angles, htAr, bfAr, volSize, Res, CropHalfWidth, cylinderRad, cylinderAtomBufferRad, fixedfa_small, projections, calcR_halfRad, useInd, SF)
element_flags = [1, 2];
R_elements = zeros(1, numel(element_flags));
parfor element_flag = element_flags
    R_elements(element_flag) = getRByElement_Cylinder(model, curr_atoms_INC, curr_IncInd_i, element_flag, angles, htAr, bfAr, volSize, Res, CropHalfWidth, cylinderRad, cylinderAtomBufferRad, fixedfa_small, projections, calcR_halfRad, useInd, SF)
end

end