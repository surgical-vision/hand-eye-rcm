function [XL, XR] = HandEye(camera_pose, robot_pose, FLAG_ALGO, stereoT, RCM_robot, RCM_grid, points, K)
    N = size(camera_pose, 3);
    
    h = optimoptions(@lsqnonlin, 'Algorithm', 'levenberg-marquardt', ...
    'MaxFunEvals', 10000, 'Display', 'Off');
    
    camRTcamL = stereoT;
    
    AR = zeros(4, 4, 4*N);
    AL = zeros(4, 4, 4*N);
    B = zeros(4, 4, 4*N);
    for i = 1:N
        if (i == N)
            AL(:, :, i) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 1);
            AR(:, :, i) = camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 2);
            AL(:, :, i + N) = camRTcamL\camera_pose(:, :, 1, 2)/(camera_pose(:, :, i, 2))*camRTcamL;
            AR(:, :, i + N) = camRTcamL*camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 1)/camRTcamL;
            AL(:, :, i + 2*N) = camRTcamL\camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 1);
            AR(:, :, i + 2*N) = camRTcamL*camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 2);
            AL(:, :, i + 3*N) = camera_pose(:, :, 1, 1)/camera_pose(:, :, i, 2)*camRTcamL;
            AR(:, :, i + 3*N) = camera_pose(:, :, 1, 2)/camera_pose(:, :, i, 1)/camRTcamL;
            B(:, :, i) = robot_pose(:, :, 1)\robot_pose(:, :, i);
        else
            AL(:, :, i) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 1);
            AR(:, :, i) = camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 2);
            AL(:, :, i + N) = camRTcamL\camera_pose(:, :, i + 1, 2)/(camera_pose(:, :, i, 2))*camRTcamL;
            AR(:, :, i + N) = camRTcamL*camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 1)/camRTcamL;
            AL(:, :, i + 2*N) = camRTcamL\camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 1);
            AR(:, :, i + 2*N) = camRTcamL*camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 2);
            AL(:, :, i + 3*N) = camera_pose(:, :, i + 1, 1)/camera_pose(:, :, i, 2)*camRTcamL;
            AR(:, :, i + 3*N) = camera_pose(:, :, i + 1, 2)/camera_pose(:, :, i, 1)/camRTcamL;
            B(:, :, i) = robot_pose(:, :, i + 1)\robot_pose(:, :, i);
        end
    end
    B(:, :, N + 1:2*N) = B(:, :, 1:N);
    B(:, :, 2*N + 1:3*N) = B(:, :, 1:N);
    B(:, :, 3*N + 1:4*N) = B(:, :, 1:N);

    switch FLAG_ALGO
        case 1 %% RCM algorithm %%
            XR = HandEye_RCM(AR, B, RCM_grid, RCM_robot(1:3, 4), camera_pose(:, :, :, 2), robot_pose);
            XL = HandEye_RCM(AL, B, RCM_grid, RCM_robot(1:3, 4), camera_pose(:, :, :, 1), robot_pose);
        case 2 %% Tsai's algorithm %%
            XR = HandEye_Tsai(AR, B);
            XL = HandEye_Tsai(AL, B);
        case 3 %% IDQ algorithm %%
            XR = HandEye_IDQ(AR, B);
            XL = HandEye_IDQ(AL, B);
        case 4 %% reproj algorithm %%
            XR = HandEye_Tsai(AR, B);
            XL = HandEye_Tsai(AL, B);
            
            Zbuffer = zeros(4, 4, N);
            
            for i = 1:N
                
                Zbuffer(:, :, i) = (camera_pose(:, :, i, 2) \ XR ) / robot_pose(:, :, i);
                
            end
            Z = avgTransformation(Zbuffer);
            XZR = [rodrigues(XR(1:3, 1:3)); XR(1:3, 4); rodrigues(Z(1:3, 1:3)); Z(1:3, 4)];
            
            for i = 1:N
                
                Zbuffer(:, :, i) = (camera_pose(:, :, i, 1) \ XL ) / robot_pose(:, :, i);
                
            end
            Z = avgTransformation(Zbuffer);
            
            XZL = [rodrigues(XL(1:3, 1:3)); XL(1:3, 4); rodrigues(Z(1:3, 1:3)); Z(1:3, 4)];
            refine_XZL = lsqnonlin(@(x) refine_reproj(x, robot_pose, points, K, 1, camera_pose(:, :, :, 2)), ...
                        XZL, [], [], h);
            refine_XZR = lsqnonlin(@(x) refine_reproj(x, robot_pose, points, K, 2, camera_pose(:, :, :, 2)), ...
                        XZR, [], [], h);
            XL = [rodrigues(refine_XZL(1:3)), refine_XZL(4:6);0 0 0 1];
            XR = [rodrigues(refine_XZR(1:3)), refine_XZR(4:6);0 0 0 1];
        case 5 %% absor algorithm %%
            
            XR = HandEye_Tsai(AR, B);
            XL = HandEye_Tsai(AL, B);
            
            refine_XR = lsqnonlin(@(x) optimRCM(x, camera_pose(:, :, :, 2), robot_pose, RCM_robot(1:3, 4), RCM_grid), ...
                        [rodrigues(XR(1:3, 1:3)); XR(1:3, 4)], [], [], h);
            refine_XL = lsqnonlin(@(x) optimRCM(x, camera_pose(:, :, :, 1), robot_pose, RCM_robot(1:3, 4), RCM_grid), ...
                        [rodrigues(XL(1:3, 1:3)); XL(1:3, 4)], [], [], h);
            XL = [rodrigues(refine_XL(1:3)), refine_XL(4:6);0 0 0 1];
            XR = [rodrigues(refine_XR(1:3)), refine_XR(4:6);0 0 0 1];
            
        otherwise
            error('Invalid ALGORITHM FLAG.');
    end
    

%     
%     if (FLAG_ALGO == 2) %%Optimise using RCM info#
%         xr_init = [rodrigues(XR(1:3, 1:3)); XR(1:3, 4)];
%         xl_init = [rodrigues(XL(1:3, 1:3)); XL(1:3, 4)];
%         [refine_XR, ~, ~, ~, ~, ~, ~] = lsqnonlin(@(xr) optimRCM(xr, camera_pose(:, :, :, 2), robot_pose, RCM_robot(1:3, 4), RCM_grid),...
%                                         xr_init, [], [], h);
%         [refine_XL, ~, ~, ~, ~, ~, ~] = lsqnonlin(@(xl) optimRCM(xl, camera_pose(:, :, :, 1), robot_pose, RCM_robot(1:3, 4), RCM_grid),...
%                                         xl_init, [], [], h);
%                                     
%         XR = [rodrigues(refine_XR(1:3)), refine_XR(4:6);0 0 0 1];
%         XL = [rodrigues(refine_XL(1:3)), refine_XL(4:6);0 0 0 1];
%     elseif (FLAG_ALGO == 3) || (FLAG_ALGO == 6) 
%         xr_init = [rodrigues(XR(1:3, 1:3)); XR(1:3, 4)];
%         xl_init = [rodrigues(XL(1:3, 1:3)); XL(1:3, 4)];
%         [refine_XR, ~, ~, ~, ~, ~, ~] = lsqnonlin(@(xr) optimAXXB(xr, AR, B),...
%                                         xr_init, [], [], h);
%         [refine_XL, ~, ~, ~, ~, ~, ~] = lsqnonlin(@(xl) optimAXXB(xl, AL, B),...
%                                         xl_init, [], [], h);
%                                     
%         XR = [rodrigues(refine_XR(1:3)), refine_XR(4:6);0 0 0 1];
%         XL = [rodrigues(refine_XL(1:3)), refine_XL(4:6);0 0 0 1];
%     end
    
end