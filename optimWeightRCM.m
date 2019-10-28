function f = optimWeightRCM(w, K, K_extend)

    newK = [(1 - w)*K; (w)*K_extend];
    [~, ~, V] = svd(newK);
    f = newK*V(:, 4);
    
end