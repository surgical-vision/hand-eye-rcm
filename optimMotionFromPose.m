function [A, B] = optimMotionFromPose(camera_pose, robot_pose, FLAG_STEREO, stereoT)
    N = size(robot_pose, 3);
    pair_selection = zeros(N, 2);
    pair_selection(:, 1) = (1:N)';
    
    for i = 1:N
        buf = zeros(1, N);
        buf(i) = -999;
        for j = setdiff(1:N, i)
            r1 = rodrigues(camera_pose(1:3, 1:3, i, 1));
            r2 = rodrigues(camera_pose(1:3, 1:3, j, 1));
            th1 = acos(r1'*r2/(norm(r1)*norm(r2)));
            relative_rot = camera_pose(1:3, 1:3, j, 1)/camera_pose(1:3, 1:3, i, 1);
            th2 = norm(rodrigues(relative_rot));
            buf(j) = sin(th1) + th2;
        end
        [~, ind] = max(buf);
        pair_selection(i, 2) = ind;
    end
    
    if (FLAG_STEREO)
        A = zeros(4, 4, 2*N);
        B = zeros(4, 4, 2*N);
        
        for i = 1:N
            A(:, :, i) = camera_pose(:, :, pair_selection(i, 2), 1)/camera_pose(:, :, i, 1);
            A(:, :, i + N) = stereoT*camera_pose(:, :, pair_selection(i, 2), 2)/camera_pose(:, :, i, 2)/stereoT;
            B(:, :, i) = robot_pose(:, :, pair_selection(i, 2))\robot_pose(:, :, i);
        end
        
        B(:, :, i+1:2*N) = B(:, :, 1:N);
    else
        
        A = zeros(4, 4, N);
        B = zeros(4, 4, N);
        
        for i = 1:N
            A(:, :, i) = camera_pose(:, :, pair_selection(i, 2), 1)/camera_pose(:, :, i, 1);
            B(:, :, i) = robot_pose(:, :, pair_selection(i, 2))\robot_pose(:, :, i);
        end
        
    end
end