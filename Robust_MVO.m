function  x_optimal = Robust_MVO(mu, Q, db_size)

    Q_gb = [Q,zeros(size(mu, 1),1);zeros(1,size(mu, 1) + 1)];
    Theta = diag(diag(Q))/db_size;
    epsilon = sqrt(chi2inv(0.9, size(mu, 1)));
    A = [ones(1, size(mu, 1)), 0 ]; 
    obj = [-mu', epsilon];
    b = ones(1, 1);
    Qua = [Theta,zeros(size(mu, 1),1);zeros(1,size(mu, 1) + 1)];
    Qua(size(mu, 1) + 1, size(mu, 1) + 1) = -1;
    v_type = blanks(size(mu, 1) + 1);
    v_type(1:end) = "C";
    lb = -10000 * ones(1, size(mu, 1) + 1);
    lb(1,size(mu, 1) + 1) = 0;

    clear model;
    model.Q = sparse(50 * Q_gb);
    model.A = sparse(A);
    model.obj = obj;
    model.sense = '=';
    model.rhs = b;
    model.vtype = v_type;
    model.lb = lb;

    model.quadcon(1).Qc = sparse(Qua);
    model.quadcon(1).q = zeros(size(mu, 1) + 1,1);
    model.quadcon(1).rhs = 0;
    model.quadcon(1).sense = '=';
    
    result = gurobi(model);

    x_optimal = result.x(1:size(mu, 1), :);

end