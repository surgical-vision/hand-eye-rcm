function err = refine_reproj(x, baseTrobot, points, K, LRindex, camTgridt)

    camTrobot = [rodrigues(x(1:3)), x(4:6);0 0 0 1];
    gridTbase = [rodrigues(x(7:9)), x(10:12);0 0 0 1];
    
    if (LRindex == 1)
        KK.fc = K.fcL;
        KK.cc = K.ccL;
        KK.kcbp = K.kcbpL;    
    elseif (LRindex == 2)
        KK.fc = K.fcR;
        KK.cc = K.ccR;
        KK.kcbp = K.kcbpR;
    end
    
    N = size(baseTrobot, 3);
    
    err_val = zeros(N, 1);
    
    for i = 1:N
        
        camTgrid = camTrobot / (gridTbase * baseTrobot(:, :, i));
        [xd, ~, ~, ~, ~, ~, ~] = project_points2(points.gt(1:3, :), rodrigues(camTgrid(1:3, 1:3)), camTgrid(1:3, 4), KK.fc, KK.cc, KK.kcbp, 0);
        [xdp, ~, ~, ~, ~, ~, ~] = project_points2(points.gt(1:3, :), rodrigues(camTgridt(1:3, 1:3, i)), camTgridt(1:3, 4, i), KK.fc, KK.cc, KK.kcbp, 0);
        err = xdp - xd;
        
        err_val(i) = mean(sqrt(sum(err.^2, 1)));
    end
    
    err = sqrt(mean(err));
end