%% MIE377 (Winter 2018) - Project 2
%
% Student Name: Shoujun Feng
% Student ID:   1001117506

clc
clear all
format long

% Program Start
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 1. Read input files
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the stock weekly prices and factors weekly returns
adjClose = readtable('Project2_Data_adjClose.csv', 'ReadRowNames', true);
adjClose.Properties.RowNames = cellstr( ...
                            datetime(adjClose.Properties.RowNames, ...
                            'InputFormat','uuuu/MM/dd'));

factorRet = readtable('Project2_Data_FF_factors.csv', 'ReadRowNames', true);
factorRet.Properties.RowNames = cellstr( ...
                            datetime(factorRet.Properties.RowNames, ...
                            'InputFormat','uuuu/MM/dd'));

riskFree = factorRet(:,4);
factorRet = factorRet(:,1:3);

% Identify the tickers and the dates 
tickers = adjClose.Properties.VariableNames';
dates   = datetime(factorRet.Properties.RowNames);

% Calculate the stocks' weekly EXCESS returns
prices  = table2array(adjClose);
returns = ( prices(2:end,:) - prices(1:end-1,:) ) ./ prices(1:end-1,:);
returns = returns - ( diag( table2array(riskFree) ) * ones( size(returns) ) );
returns = array2table(returns);
returns.Properties.VariableNames = tickers;
returns.Properties.RowNames = cellstr(datetime(factorRet.Properties.RowNames));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2. Define your initial parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of shares outstanding per company
mktNoShares = [36.40 9.86 15.05 10.91 45.44 15.43 43.45 33.26 51.57 45.25 ...
                69.38 8.58 64.56 30.80 36.92 7.70 3.32 54.68 38.99 5.78]' * 1e8;

% Initial budget to invest
initialVal = 100 * ones(1, 5);

% Start of in-sample calibration period 
calStart = datetime('2012-01-01');
calEnd   = calStart + calmonths(12) - days(1);

% Start of out-of-sample test period 
testStart = datetime('2013-01-01');
testEnd   = testStart + calmonths(6) - days(1);

% Number of investment periods (each investment period is 6 months long)
NoPeriods = 6;

% Investment strategies
% Note: You must populate the functios MVO.m, MVO_card.m and BL.m with your
% own code to construct the optimal portfolios. 
funNames  = {'MVO' 'Robust MVO confidence = 0.9' 'Re sampling' 'most deverse' 'CVaR'};
%funNames  = {'MVO' 'MVO (Card=12)' 'B-L Model'};
NoMethods = length(funNames);

funList = {'MVO' 'Robust_MVO' 'Re_sampling' 'most_diverse' 'CVaR'};
% funList = {'MVO' 'MVO_card' 'BL'};
funList = cellfun(@str2func, funList, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 3. Construct and rebalance your portfolios
%
% Here you will estimate your input parameters (exp. returns, cov. matrix,
% etc) from the Fama-French factor models. You will have to re-estimate 
% your parameters at the start of each rebalance period, and then 
% re-optimize and rebalance your portfolios accordingly. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initiate counter for the number of observations per investment period
toDay = 0;

for t = 1 : NoPeriods
    
    % Subset the returns and factor returns corresponding to the current
    % calibration period.
    disp('period');
    disp(t);
    disp('calStart:');
    disp(calStart);
    disp('calEnd:');
    disp(calEnd);
    disp('testStart');
    disp(testStart);
    periodReturns = table2array( returns( calStart <= dates & dates <= calEnd, :) );
    periodFactRet = table2array( factorRet( calStart <= dates & dates <= calEnd, :) );
    currentPrices = table2array( adjClose( ( calEnd - days(7) ) <= dates ... 
                                                    & dates <= calEnd, :) )';
    
    % Subset the prices corresponding to the current out-of-sample test 
    % period.
    periodPrices = table2array( adjClose( testStart <= dates & dates <= testEnd,:) );
    
    % Set the initial value of the portfolio or update the portfolio value
    if t == 1
        
        currentVal(t,:) = initialVal;
        
    else
        for i = 1 : NoMethods
            
            currentVal(t,i) = currentPrices.' * NoShares{i};
            
        end
    end
    
    % Update counter for the number of observations per investment period
    fromDay = toDay + 1;
    toDay   = toDay + size(periodPrices,1);
    
    % *************** WRITE YOUR CODE HERE ***************
    %----------------------------------------------------------------------
    
    % Calculate your initial exp. returns and covariance matrix using the
    % Fama-French factor model. You may find these values as you prefer,
    % the code given below is just an example. 
    
    % B is the period factor matrix used in OLS correlation
    B = [ones(size(periodFactRet, 1), 1), periodFactRet];
    
    % Use the OLS closed form solution to get coefficients
    coeffi = (B.' * B) \ B.' * periodReturns;
    
    % alphas for all stocks
    A = coeffi(1, :).';          % n x 1 vector of alphas
    
    % betas for all stocks
    V = coeffi(2:end, :);          % m x n matrix of betas
    
    f_bar = geomean(1 + periodFactRet) - 1;      % m x 1 vector of factor expected returns
    F = cov(periodFactRet);                        % m x m factor covariance matrix
    
    epsilon = periodReturns - B * coeffi;   % Regression residuals
    D = diag(var(epsilon));          % Diagonal n x n matrix of residual variance
    
    mu = A + V' * f_bar';         % n x 1 vector of asset exp. returns
    Q = V' * F * V + D;        % n x n asset covariance matrix
    
    %----------------------------------------------------------------------
    
    % Define the target return for the 2 MVO portfolios
    targetRet = mean(mu);
    
    % Optimize your portfolios to get the weights 'x'
    % Note: You need to write the code for the 3 functions ('MVO.m', 
    % 'MVO_card.m' and 'BL.m') to find your optimal portfolios
    
    x{1}(:,t) = funList{1}(mu, Q, targetRet); 
    x{2}(:,t) = funList{2}(mu, Q, size(periodReturns, 1)); 
    x{3}(:,t) = funList{3}(mu, Q, targetRet);
    x{4}(:,t) = funList{4}(mu, Q, targetRet, 12);
    x{5}(:,t) = funList{5}(mu, Q, targetRet, currentPrices);
    
    % Calculate the optimal number of shares of each stock you should hold
    for i = 1 : NoMethods
        
        % Number of shares your portfolio holds per stock
        NoShares{i} = x{i}(:,t) .* currentVal(t,i) ./ currentPrices;
        
        % Weekly portfolio value during the out-of-sample window
        portfValue(fromDay:toDay,i) = periodPrices * NoShares{i};
        
        % *************** WRITE YOUR CODE HERE ***************
        %------------------------------------------------------------------
        
        % Calculate your transaction costs for the current rebalance
        % period. The first period does not have any cost since you are
        % constructing the portfolios for the 1st time. 
        
        % Calculate the transaction cost
        if t ~= 1
           
            tCost(t-1, i) = 0.005 * currentPrices' *  NoShares{i};
            
        end
        
        NoSharesOld{i} = NoShares{i};
        %------------------------------------------------------------------
        
    end

    % Update your calibration and out-of-sample test periods
    calStart = calStart + calmonths(6);
    calEnd   = calStart + calmonths(12) - days(1);
    
    testStart = testStart + calmonths(6);
    testEnd   = testStart + calmonths(6) - days(1);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 4. Results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% *************** WRITE YOUR CODE HERE ***************
%--------------------------------------------------------------------------

% Calculate the portfolio average return, variance (or standard deviation),
% or any other performance and/or risk metric you wish to include in your
% report.

% Calculate portfolio weekly return
portfolio_return = portfValue(2:end, :) ./ portfValue(1:end - 1, :) - 1;

% Calculate the average of the weekly return by geo_mean
portf_avg_return = geo_mean(1 + portfolio_return) - 1;

% Calculate the portfolio variance
portf_variance = var(portfolio_return);

% Calculate the total transaction cost for portfolios
total_transaction_cost = sum(tCost);

% Display portfolio results
disp('Portfolio result summary')
disp('Average portfolio returns are:');
disp(portf_avg_return);
disp('Portfolio variance:');
disp(portf_variance);
disp('Portfolio initial value:');
disp(initialVal);
disp('Portfolio final value:');
disp(portfValue(end, :));
disp('Total transaction cost:');
disp(total_transaction_cost );
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% 4.1 Plot the portfolio values 
%--------------------------------------------------------------------------
plotDates = dates(datetime('2013-01-01') <= dates);

fig1 = figure(1);
plot(plotDates, portfValue(:,1))
hold on
plot(plotDates, portfValue(:,2))
hold on
plot(plotDates, portfValue(:,3))
hold on
plot(plotDates, portfValue(:,4))
hold on
plot(plotDates, portfValue(:,5))
legend(funNames, 'Location', 'eastoutside','FontSize',12);
datetick('x','dd-mmm-yyyy','keepticks','keeplimits');
set(gca,'XTickLabelRotation',30);
title('Portfolio Values of three models', 'FontSize', 14)
ylabel('Value','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig1,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig1,'Position');
set(fig1,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig1,'fileName','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig1,'Portfolio values comparison','-dpng','-r0');

%--------------------------------------------------------------------------
% 4.2 Plot the portfolio weights 
%--------------------------------------------------------------------------

% MVO Plot
fig2 = figure(2);
area(x{1}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights of MVO', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig2,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig2,'Position');
set(fig2,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig2,'fileName2','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig2,'MVO_weight','-dpng','-r0');


% Robust MVO
fig3 = figure(3);
area(x{2}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights of Robust MVO', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig3,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig3,'Position');
set(fig3,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig3,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig3,'Robust_MVO_weight','-dpng','-r0');


% Resampling Plot
fig4 = figure(4);
area(x{3}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights of Resampling', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig4,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig4,'Position');
set(fig4,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig4,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig4,'Re_sampling_weight','-dpng','-r0');

% Most_diverse Plot
fig5 = figure(5);
area(x{4}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights of Most diverse', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig5,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig5,'Position');
set(fig5,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig4,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig5,'Most_diverse_weight','-dpng','-r0');

% CVaR Plot
fig6 = figure(6);
area(x{5}')
legend(tickers, 'Location', 'eastoutside','FontSize',12);
title('Portfolio Weights of CVaR', 'FontSize', 14)
ylabel('Weights','interpreter','latex','FontSize',12);
xlabel('Rebalance Period','interpreter','latex','FontSize',12);

% Define the plot size in inches
set(fig6,'Units','Inches', 'Position', [0 0 8, 5]);
pos1 = get(fig6,'Position');
set(fig6,'PaperPositionMode','Auto','PaperUnits','Inches',...
    'PaperSize',[pos1(3), pos1(4)]);

% If you want to save the figure as .pdf for use in LaTeX
% print(fig4,'fileName3','-dpdf','-r0');

% If you want to save the figure as .png for use in MS Word
print(fig6,'C_VaR_weight','-dpng','-r0');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Program End