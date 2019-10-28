function f = optimAXXBw(x, A, B, w)
    N = size(A, 3);
    X = [rodrigues(x(1:3)), x(4:6);0 0 0 1];
    E = zeros(6, N);  
    for i = 1:N
        LHS = logm(A(:, :, i)*X);
        RHS = logm(X*B(:, :, i));
        d = LHS - RHS;
        E(:, i) = abs([d(3, 2); d(3, 1); d(2, 1); d(1:3, 4)]);
    end
    E = sum(E, 2);
    f = E'*diag(w)*E;
end