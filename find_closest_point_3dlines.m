function f = find_closest_point_3dlines(x, origin, d)

    N = size(origin, 2);
    
    f = 0;
    tt = zeros(1, N);
    for i = 1:N
        proj_x = d(:, i)'*(x - origin(:, i));
        if (proj_x > 0)
            vec = x - (origin(:, i) + proj_x*d(:, i));
            f = f + vec'*vec;
            tt(i) = vec'*vec;
        else
            vec = x - origin(:, i);
            f = f + vec'*vec;
            tt(i) = vec'*vec;
        end
    end

end