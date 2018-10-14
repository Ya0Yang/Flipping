% addpath ../../
tic
rng(42)

% Modified H and B from Spatial realignment

htAr = [0.4140 1.5446 0];
bfAr = [5.2121 6.3684 0.1];
           
SF_Iguess = 1;

curr_measurementnum=1;
SPHyn = 1;
%%
inputDir = '';

saveprefix = 'output';

pj_filename              = [inputDir 'multislice_NiPt_171118_tomo2_BA_gauss_NegDef_shifted_norm_sameAsLinear.mat'];
% pj_filename              = [inputDir 'multislice_NiPt_171118_tomo2_BA_gauss_NegDef_shifted_09_norm_sameAsLinear.mat'];
% pj_filename              = [inputDir 'Linear_NiPt_1118t2_Yao_235_1010.mat'];

angle_filename           = [inputDir 'Angles_Multislice_NiPt_171118_tomo2_BA_gauss_NegDef.mat'];
atomtype_filename        = [inputDir 'NiPt_1118t2_model_atomtype.mat'];
model_filename           = [inputDir 'NiPt_1118t2_model_OriOri_1009.mat'];
IncInd_filename          = [inputDir 'NiPt_1118t2_shell_ind_rand1000_1009.mat'];
% IncInd_filename          = [inputDir 'NiPt_1118t2_shell_ind_rand300_1012.mat'];
% IncInd_filename          = [inputDir 'Linear_rand_atomType_NiPt_171118_tomo2.mat'];

% IncInd_filename          = [inputDir 'NiPt_1118t2_inside_ind_0920.mat'];
%%

Allangles = importdata(angle_filename);
angles = Allangles(:,:);

model = importdata(model_filename);

atoms = importdata(atomtype_filename);

projections = importdata(pj_filename);

IncInd      = importdata(IncInd_filename);
% IncInd_arr = randperm(length(atoms));
% IncInd     = IncInd_arr(1:500);

%%
projections = My_paddzero(projections,[size(projections,1)+20 size(projections,2)+20 size(projections,3)]);

% projections = projections(2:end,2:end,:);
%%
atoms1 = atoms;
for i = 1:round(length(IncInd)*2/3)
    rand_type = randperm(3);
    atoms1(IncInd(i)) = rand_type(1);
end
% save(sprintf('%s_ori_rand_changed_atomtype.mat',saveprefix),'atoms','atoms1','IncInd');

%%
atoms_INC = atoms1;

% volSize= 234;
volSize= size(projections,1);


Res = 0.469;
CropHalfWidth = 4; % in Pixel
cylinderRad = 4; % in Angstrom

cylinderAtomBufferRad = 1; % in Angstrom

calcR_halfRad = 5;

fixedfa_big = make_fixedfa(volSize, Res);
fixedfa_small = make_fixedfa(1+4*round(cylinderRad/Res), Res);

[X,Y,~] = ndgrid(-calcR_halfRad:calcR_halfRad,-calcR_halfRad:calcR_halfRad,1:size(angles,1));
if SPHyn == 1
  useInd = find( (X.^2+Y.^2) < (calcR_halfRad+0.5)^2);
else
  useInd = 1:length(X(:));
end
clear X Y

%%

Atomtype_INC = atoms_INC;

curr_atoms_INC = Atomtype_INC;
prev_atoms_INC = Atomtype_INC;

iterNum=0;
endFlag = 0;
while ~endFlag
    iterNum = iterNum+1;

    curr_Projs = My_create_volProjs_from_model_exact_HB_fixedfa_fainput(model,curr_atoms_INC,htAr,bfAr,volSize, Res, CropHalfWidth, angles, fixedfa_big);
%     curr_Projs = My_paddzero(curr_Projs,[size(curr_Projs,1)+20 size(curr_Projs,2)+20 size(curr_Projs,3)]);
    fun = @(x,xdata)sum(abs(xdata(:,1) - x*xdata(:,2)));
    x0 = SF_Iguess;
    SF = lsqcurvefit(fun,x0,[projections(:) curr_Projs(:)],0)

    localR_ar_Pt = zeros(1,length(IncInd));
    localR_ar_Ni = zeros(1,length(IncInd));
    localR_ar_NA = zeros(1,length(IncInd));
    
    curr_IncInd = IncInd(randperm(length(IncInd)));

    %%CHANGE THIS BACK
    for i=1:length(curr_IncInd)     
%     for i=1:250
        curr_atoms_INC(curr_IncInd(i)) = 1;
        [cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D_fainput_parfor(model, curr_atoms_INC, curr_IncInd(i),...
            angles, [htAr],[bfAr], volSize, Res,CropHalfWidth, round(cylinderRad/Res), round(cylinderAtomBufferRad/Res), fixedfa_small);

        subProjs = zeros(size(cylinerIndsAr,2),size(cylinerIndsAr,2),size(cylinerIndsAr,3));
        for j=1:size(angles,1)
          subProjs(:,:,j) = projections(cylinerIndsAr(1,:,j),cylinerIndsAr(2,:,j),j);
        end

        subProjs_crop = My_stripzero_ajp(subProjs ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        cylinderProjAr_crop = My_stripzero_ajp(cylinderProjAr ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        currR_Ni = calcR_nonorm_YY(subProjs_crop(useInd),cylinderProjAr_crop(useInd)*SF);
        
        curr_atoms_INC(curr_IncInd(i)) = 2;
        [cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D_fainput_parfor(model, curr_atoms_INC, curr_IncInd(i),...
            angles, [htAr],[bfAr], volSize, Res,CropHalfWidth, round(cylinderRad/Res), round(cylinderAtomBufferRad/Res), fixedfa_small);

        subProjs = zeros(size(cylinerIndsAr,2),size(cylinerIndsAr,2),size(cylinerIndsAr,3));
        for j=1:size(angles,1)
          subProjs(:,:,j) = projections(cylinerIndsAr(1,:,j),cylinerIndsAr(2,:,j),j);
        end

        subProjs_crop = My_stripzero_ajp(subProjs ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        cylinderProjAr_crop = My_stripzero_ajp(cylinderProjAr ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        currR_Pt = calcR_nonorm_YY(subProjs_crop(useInd),cylinderProjAr_crop(useInd)*SF);
        
        % Non-atom
        curr_atoms_INC(curr_IncInd(i)) =3;
        [cylinderProjAr, cylinerIndsAr] = My_create_cylinder_proj_from_model_2D_fainput_parfor(model, curr_atoms_INC, curr_IncInd(i),...
            angles, [htAr],[bfAr ],volSize, Res, CropHalfWidth, round(cylinderRad/Res), round(cylinderAtomBufferRad/Res), fixedfa_small);
        subProjs = zeros(size(cylinerIndsAr,2),size(cylinerIndsAr,2),size(cylinerIndsAr,3));
        for j=1:size(angles,1)
          subProjs(:,:,j) = projections(cylinerIndsAr(1,:,j),cylinerIndsAr(2,:,j),j);
        end

        subProjs_crop = My_stripzero_ajp(subProjs ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        cylinderProjAr_crop = My_stripzero_ajp(cylinderProjAr ,[calcR_halfRad*2+1 calcR_halfRad*2+1 size(angles,1)]);
        currR_NA = calcR_nonorm_YY(subProjs_crop(useInd),cylinderProjAr_crop(useInd)*SF);
        
        localR_ar_Ni(i) = currR_Ni;
        localR_ar_Pt(i) = currR_Pt;
        localR_ar_NA(i) = currR_NA;

%         [~,minInd] = min([currR_Ni, currR_Pt]);
        [~,minInd] = min([currR_Ni, currR_Pt, currR_NA]);
        curr_atoms_INC((curr_IncInd(i))) = minInd;
        
        if mod(i,100) == 0
            disp(i)
        end
    end
    
    localR_ar = [localR_ar_Ni;localR_ar_Pt;localR_ar_NA];
%     localR_ar = [localR_ar_Ni;localR_ar_Pt];
%     save(sprintf('%s_Rar_iter_%d.mat',saveprefix,iterNum),'localR_ar','curr_atoms_INC','SF','curr_IncInd');
    
    if sum(prev_atoms_INC~=curr_atoms_INC)==0
        fprintf(1,'iter %d, converged!\n',iterNum);
        endFlag=1;
    else
        fprintf(1,'iter %d, %d atoms changed!\n',iterNum, sum(prev_atoms_INC~=curr_atoms_INC));
        prev_atoms_INC=curr_atoms_INC;
        endFlag=0;
        
    end
end
toc

% flipped=find(atoms-curr_atoms_INC~=0)
% curr_atoms_INC((atoms-curr_atoms_INC~=0))
