function X = HandEye_RCM(A, B, RCM_grid, RCM_base, camera_pose, robot_pose)
    h = optimoptions(@fminunc, 'MaxFunEvals', 20000, 'Display', 'Off');
    h_con = optimoptions(@fmincon, 'MaxFunEvals', 10000, 'Display', 'Off');

    N = size(A, 3);
    N_Transform = size(robot_pose, 3);

    RCM_robot = zeros(4, N_Transform);
    RCM_cam = zeros(4, N_Transform);
    
    X = zeros(4, 4);
    X(4, 4) = 1;
    v_a = zeros(3, N);
    v_b = zeros(3, N);
    LHS1 = zeros(3*N, 3);

    
    for j = 1:N
        LHS1(3*j - 2:3*j, 1:3) = A(1:3, 1:3, j) - eye(3);
        a = logm(A(:, :, j));
        b = logm(B(:, :, j));
        v_a(:, j) = a(1:3, 4);
        v_b(:, j) = b(1:3, 4);
    end
    
    for i = 1:N_Transform
        RCM_robot(:, i) = robot_pose(:, :, i)\[RCM_base; 1];
        RCM_cam(:, i) = camera_pose(:, :, i)*[RCM_grid; 1];
    end
    
    [optim_rot, valR] = fmincon(@(x) optim_rot_rcm2(x, RCM_cam, RCM_robot, LHS1(:, 1:3), A, B), ...
       [0, 0, pi/180], [], [], [], [], [], [], @(x) mycon(x, robot_pose), h_con);
    
%     [optim_rot, valR] = fminunc(@(x) optim_rot_rcm2(x, RCM_cam, RCM_robot, LHS1(:, 1:3), A, B), ...
%         zeros(3, 1), h);

    X(1:3, 1:3) = rodrigues(optim_rot);
    
    RCM_diff = RCM_cam(1:3, :) - X(1:3, 1:3)*RCM_robot(1:3, :);
    RHS1 = X(1:3, 1:3)*reshape(B(1:3, 4, :), 3, N) - reshape(A(1:3, 4, :), 3, N);
    RHS1 = RHS1(:);

    [X(1:3, 4), valT] = fminunc(@(x) optim_tran_rcm(x, RCM_diff, LHS1(:, 1:3), RHS1), ...
                                    zeros(3, 1), h);

end

function [c, ceq] = mycon(x, robot_pose)

    N = size(robot_pose, 3);
    ceq = 0;
    c = zeros(N + 1, 1);
    
    R = rodrigues(x);
    
    for i = 1:N
        
        worldRcam = robot_pose(1:3, 1:3, i)*R';
        theta = real(acos(worldRcam(1:3, 3)'*robot_pose(1:3, 3, i)));
        c(i) = theta - 4*pi/180;
        
    end
    
    c(N + 1) = norm(x) - pi;
    
end

function f = optim_tran_rcm(x, diff_RCM, LHS, RHS)

    t = x;
    
    val1 = repmat(t, 1, size(diff_RCM, 2)) - diff_RCM;
    val1 = sqrt(diag(val1'*val1));
    
%     X = [R, t; 0 0 0 1];
%     
%     N = size(camera_pose, 3);
%     
%     f = sum(val1) + std_x^2 + std_y^2 + std_z^2;
    val2 = LHS*t - RHS(:);
    val2 = reshape(val2, 3, numel(val2)/3);
    val2 = sqrt(diag(val2'*val2));

    f = sum(val1) + sum(val2);
end


function f = optim_rot_rcm2(x, RCM_cam, RCM_robot, LHS, A, B)

    R = rodrigues(x);
    
    N = size(A, 3);

    t_RCM = mean(RCM_cam(1:3, :) - R*RCM_robot(1:3, :), 2);

    RHS = R * reshape(B(1:3, 4, :), 3, N) - reshape(A(1:3, 4, :), 3, N);
    t_HE = LHS\RHS(:);
    
    cen_RCM_cam = RCM_cam(1:3, :) - repmat(mean(RCM_cam(1:3, :), 2), 1, N/4);
    cen_RCM_robot = RCM_robot(1:3, :) - repmat(mean(RCM_robot(1:3, :), 2), 1, N/4);
    
    err_absor = cen_RCM_cam - R * cen_RCM_robot;
    
    f = norm(t_RCM - t_HE) + sum(vecnorm(err_absor));
end