function support_index = OMP(yt,dict,num_iteration)
    residual = yt;
    support_index = [];
    
    for iter = 1:num_iteration
        corr = dict' * residual;
        [~, idx] = max(abs(corr));
        support_index = unique([support_index, idx]);
        A_sel = dict(:, support_index);
        x_hat = pinv(A_sel) * yt;
        residual = yt - A_sel * x_hat;
    end
end