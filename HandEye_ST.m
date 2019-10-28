function X = HandEye_ST(A, B)
    h2 = optimset('GradObj', 'on', 'Algorithm', 'interior-point', 'Display', 'off');

    converge_thres = 10;
    counter_rotation = 0;
    counter_translation = 0;
    max_ite = 50;
    
    N = size(A, 3);

    v_a = zeros(3, N);
    v_b = zeros(3, N);
    om_a = zeros(3, N);
    om_b = zeros(3, N);
    K = zeros(8*N, 4);

    X = HandEye_IDQ(A, B);
    
    LHS = zeros(3*N, 3);
    RHS = zeros(3*N, 1);
    

    %Construct L and Lp matrix
    L = zeros(4*N, 4);
    Lp = zeros(4*N, 4);
    
    for j = 1:N
        a = logm(A(:, :, j));
        b = logm(B(:, :, j));
        v_a(:, j) = a(1:3 ,4);
        v_b(:, j) = b(1:3, 4);
        om_a(:, j) = rodrigues(A(1:3, 1:3, j));
        om_b(:, j) = rodrigues(B(1:3, 1:3, j));
%         LHS(3*j - 2:3*j, :) = skew3(om_a(:, j));
%         LHS(3*j - 2:3*j, :) = skew3(X(1:3, 1:3)*om_b(:, j));
        [a, ap] = getDualQ(A(1:3, 1:3, j), A(1:3, 4, j));
        [b, bp] = getDualQ(B(1:3, 1:3, j), B(1:3, 4, j));
        x = a - b;
        y = a + b;
        z = ap - bp;
        w = ap + bp;
        L(4*j-3:4*j, :) = [x(4), -x(1:3)';x(1:3), skew3(y(1:3))+x(4)*eye(3)];
        Lp(4*j-3:4*j, :) = [z(4), -z(1:3)';z(1:3), skew3(w(1:3))+z(4)*eye(3)];
%         K(14*j - 13:14*j - 8, :) = [x(1:3), skew3(y(1:3)), zeros(3, 4);
%                              z(1:3), skew3(w(1:3)), x(1:3), skew3(y(1:3))];  
        K(8*j - 7:8*j - 4, 1:4) = [x(4), -(x(1:3))'; x(1:3), skew3(y) + x(4)*eye(3)];
%         K(12*j - 7:12*j - 4, 1:4) = [om_a(:, j) - om_b(:, j), skew3(om_a(:, j) + om_b(:, j)); 0, (om_b(:, j) - om_a(:, j))'];
    end

    counter = 0;
    
    while (((counter_rotation < converge_thres) || (counter_translation < converge_thres)) && (counter < max_ite))
        
        R_init = X(1:3, 1:3);
        t_init = X(1:3, 4);
        
        for j = 1:N
            vec_buf = v_a(:, j) - cross(t_init, R_init*om_b(:, j));
%             vec_buf = v_a(:, j) - cross(t_init, om_a(:, j));
            K(8*j - 3:8*j, 1) = [vec_buf - v_b(:, j); 0];
            K(8*j - 3:8*j, 2:4) = [skew3(vec_buf + v_b(:, j)); (-vec_buf + v_b(:, j))'];
        end
        
        [~, ~, v_basis] = svd(K);
        v_basis = v_basis(:, 4);

        qR = [v_basis(2:4); v_basis(1)];
        
        X(1:3, 1:3) = q2dcm(qR)';

        for j = 1:N
            LHS(3*j - 2:3*j, :) = skew3(X(1:3, 1:3)*om_b(:, j));
%             LHS(3*j - 2:3*j, :) = skew3(om_a(:, j));
            RHS(3*j - 2:3*j, 1) = X(1:3, 1:3)*v_b(:, j) - v_a(:, j);
        end


        X = [X(1:3, 1:3) LHS\RHS;0 0 0 1];
        
        diff = X/[R_init, t_init;0 0 0 1];
        
        if (norm(rodrigues(diff(1:3, 1:3))) < 1e-3)
            counter_rotation = counter_rotation + 1;
        else
            counter_rotation = 0;
        end
        
        if (norm(diff(1:3, 4)) < 1e-3)
            counter_translation = counter_translation + 1;
        else
            counter_translation = 0;
        end
        
        counter = counter + 1;
        ccc = rodrigues(X(1:3, 1:3));
        C(1:3, counter) = ccc;
        C(4:6, counter) = X(1:3, 4);
    end
        
%     LL = zeros(6, N);
%     RR = zeros(6, N);
%     
%     for i = 1:N
%         LL(:, i) = [om_a(:, i); v_a(:, i)];
%         RR(:, i) = [om_b(:, i); v_b(:, i)];
%     end
%     
%     x = fminunc(@(x) optimHandeye(x, LL, RR), [rodrigues(R_init); t_init], h2);
%     
%     X = [rodrigues(x(1:3)), x(4:6);0 0 0 1];
end