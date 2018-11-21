function [mvo]=ef(mu,Q)
%Produces the MVO efficient frontier for an expected return vector mu and
%covariance matrix Q. Each row of mvo.x is a portfolio corresponding to the
%expected return and risk in mvo.exp_ret and mvo.risk respectively.
%Ensure mu is a row vector.
%Assumes no shortselling and no upper or lower bound constraints.

options = optimoptions('quadprog','Display','off');
%options_lp = optimoptions('linprog','Display','off');

n=size(Q,1);
A = -mu;

%no shortselling
lb=zeros(n,1);

%constrain weights to sum to 1
Aeq=ones(1,n);
beq=1;

% Find minimum risk portfolio
% solve, pass empty matrices for f, A, b, and ub
x=quadprog(Q,[],[],[],Aeq,beq,lb,[],[],options);

% Minimum risk portfolio return
minret=mu*x;

% Maximum return portfolio
% solve min - portfolio return, s.t. budget & short selling
x=linprog(A,[],[],Aeq,beq,lb,[]);

% Maximum expected return
maxret= -A*x;

% for an efficient frontier of 50 points
stepsize=(maxret-minret)/49;

i=1;
for ret_targ=minret:stepsize:maxret
    b=-ret_targ;
    mvo.x(:,i)=quadprog(Q,[],A,b,Aeq,beq,lb,[],[],options);
    mvo.exp_ret(i,1)=(mu*mvo.x(:,i));
    mvo.risk(i,1)=((mvo.x(:,i)'*Q*mvo.x(:,i)))^.5;
    i=i+1;
end