function  x_optimal = Re_sampling(mu, Q, targetRet)

    L = 100;
    S = 10000;
    effi_num_pt = 50;
    num_asset = size(mu, 1);
    x_sum = zeros(num_asset, effi_num_pt);
    mu_sum = zeros(effi_num_pt, 1);

    for i = 1:L
        random_return = mvnrnd(mu, Q, S);
        sample_mean = mean(random_return);
        sample_covariance = cov(random_return);
        mvo_i = ef1(sample_mean, sample_covariance);
        x_sum = x_sum + mvo_i.x;
        mu_sum = mu_sum + mvo_i.exp_ret;
    end

    x_ave = x_sum / L;
    mu_ave = mu_sum / L;
    residual = mu_ave - targetRet * ones(effi_num_pt, 1);
    index = effi_num_pt;
    for i = 1:effi_num_pt
        if residual(i, 1) > 0
            index = i;
            break
        end
    end
    x_optimal = x_ave(:, index);
end