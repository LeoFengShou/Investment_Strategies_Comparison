function x_optimal = most_diverse( mu, Q, targetRet, card )
    n = size(Q, 1);
    rho = corrcov(Q);   
    f = zeros(n * (n + 1), 1);
    f(1 : n * n, 1) = reshape(rho, [n * n, 1]);
    intcon = [1: n * (n + 1)];
    Aeq = zeros(1 + n, n * (n + 1));
    beq = ones(1 + n, 1);
    Aeq (1, n * n + 1 : n * (n + 1)) = 1;
    beq(1, 1) = card;
    for i = 1:n
       Aeq(i + 1, n * (i - 1) + 1 : n * i) = 1; 
    end
    temp = []
    for i = 1:n
        temp = [temp; -1 * eye(n)];
    end
    A = [eye(n * n), temp];
    b = zeros(n * n , 1);
    lb = zeros(n * (n + 1), 1);
    ub = ones(n * (n + 1), 1);
    x_inte = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub);
    selected_index = [];
    for i = n * n + 1: n * (n + 1)
        if x_inte(i, 1) == 1
            selected_index = [selected_index, i - n * n];
        end
    end
    reduced_Q = zeros(size(selected_index, 2), size(selected_index, 2));
    for i = 1:size(selected_index, 2)
        for j = 1:size(selected_index, 2)
            reduced_Q(i, j) = Q(selected_index(i), selected_index(j));
        end
    end
    x_reduced_mvo = MVO(mu(selected_index, 1), reduced_Q, targetRet)
    x_optimal = zeros(n, 1);
    x_optimal(selected_index, 1) = x_reduced_mvo;
end

