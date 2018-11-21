function x_optimal = CVaR(mu, Q, targetRet, currentPrices)
    num_asset = size(mu, 1);
    rho = corrcov(Q);
    nPaths = 2000;
    L = chol(rho, 'lower');
    T = 26;
    N = 1;
    dt = T/N;
    confidence_level = 0.95;
   
    S = zeros(num_asset, N+1, nPaths);  % Matrix of simulated price paths
    % S(:, 1, :) = 100;           % Set the initial price of the asset
    for i = 1:num_asset
        for j = 1:nPaths
            S(i, 1, j) = currentPrices(i);
        end
    end
    % Generate paths
    for i = 1:nPaths
        for j = 1:N
            xi = L * randn(num_asset, 1); 
            for k = 1:num_asset
                S(k, j+1, i) = S(k, j, i) * exp( ( mu(k, 1) - 0.5 * Q (k, k) ) * dt ...
                                + sqrt(Q(k, k)) * sqrt(dt) * xi (k) );
            end
        end
    end 
    returns_sample = zeros(num_asset, nPaths); % returns_sample 20 * 2000
    for i = 1:nPaths
       for j = 1:num_asset
            returns_sample(j, i) = S(j, end, i) / S(j, 1, i) - 1;
       end
    end
    % have 2000 + 20 + 1 variables and 2000 + 2000 + 1  inequality constraints
    % and 1 equality constraint
    % 2000 z_s first, x_i second and gamma last
    f = [1 / ((1 - confidence_level) * nPaths) * ones(nPaths, 1); zeros(num_asset, 1); 1];
    A = zeros(2 * nPaths, nPaths + num_asset + 1);
    A(1:nPaths, 1:nPaths) = -1 * eye(nPaths);
    A(nPaths + 1: 2 * nPaths, 1: nPaths) = -1 * eye(nPaths);
    A(nPaths + 1: 2 * nPaths, end) = -1;
    for i = nPaths + 1:2 * nPaths
        A(i, nPaths + 1: nPaths + num_asset) = - (returns_sample(:, i - nPaths))';
    end
    Aeq = zeros(1, nPaths + num_asset + 1);
    Aeq(1, nPaths + 1: nPaths + num_asset) = 1;
    beq = ones(1, 1);
    b = zeros(2 * nPaths, 1);
    b(2 * nPaths + 1, 1) = -targetRet;
    x = linprog(f, A, b, Aeq, beq, [-1000000 * ones(2000,1); zeros(20, 1); -1000000], []);
    x_optimal = x(nPaths + 1: nPaths + num_asset, 1);
end