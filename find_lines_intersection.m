function P = find_lines_intersection(origin, direction_vector)

    N = size(origin, 2);
    
    LHS = -N*eye(3);
    RHS = zeros(3, 1);
%     RHS1 = zeros(3, 1);
%     RHS2 = zeros(3, 1);
%     RHS3 = zeros(3, 1);
    for i = 1:N
        LHS = LHS + direction_vector(:, i)*direction_vector(:, i)'; 
        RHS = RHS + (direction_vector(:, i)*direction_vector(:, i)' - eye(3)) * ...
                     origin(:, i);
%          RHS1 = RHS1 + (direction_vector(:, i)*direction_vector(:, i)' - eye(3)) * ...
%                      (origin(:, i) + direction_vector(:, i));
                 
%                  RHS2 = RHS2 + (direction_vector(:, i)*direction_vector(:, i)' - eye(3)) * ...
%                      (origin(:, i) - direction_vector(:, i));
                 
%                  RHS3 = RHS3 + (direction_vector(:, i)*direction_vector(:, i)' - eye(3)) * ...
%                      (origin(:, i) + 0.2*direction_vector(:, i));
    end
    
    P = LHS\RHS;
%     P1 = LHS\RHS1;
%     P2 = LHS\RHS2;
%     P3 = LHS\RHS3;

end