function f = optimRCM(x, camTgrid, baseTrobot, RCM_robot, RCM_grid)
    
    f = 0;

    camTrobot = [rodrigues(x(1:3)), x(4:6);0 0 0 1];
    N = size(baseTrobot, 3);
    
    for i = 1:N
        
        error_vec = camTgrid(:, :, i)*[RCM_grid; 1] - ...
            camTrobot*inv(baseTrobot(:, :, i))*[RCM_robot; 1];
        f = f + norm(error_vec);

    end

end