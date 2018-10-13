function fixedfa =  make_fixedfa(sizeX, Res)
    finalvol_summed = zeros( sizeX);

    kx = 1:size(finalvol_summed,1);
    ky = 1:size(finalvol_summed,2);

    MultF_X = 1/(length(kx)*Res);
    MultF_Y = 1/(length(ky)*Res);

    CentPos = round((size(finalvol_summed)+1)/2);
    [KX, KY] = ndgrid((kx-CentPos(1))*MultF_X,(ky-CentPos(2))*MultF_Y);
    q2 = KX.^2 + KY.^2 ;
    clear KX KY KZ

    fa78 = fatom_vector( sqrt(q2),78 );
    fa28 = fatom_vector( sqrt(q2),28 );
   
    fixedfa = 0.5*(fa78+fa28);
end