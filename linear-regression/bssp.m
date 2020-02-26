function solution = bssp(x,y,p,M,interceptFlag,predFlag,normFlag,opt_params)
%bssp Best Subset Selection Problem
%
% Perform Best Subset Selection Problem to find the minimal number of
%    predictors, p, to get linear regression of the form:
%           y = x*b + e
%       where
%           y: response/independent variables
%           x: explanatory/dependent variables
%           b: regression coefficients (to be solved for)
%           e: error/noise
% By solving a LP problem of the form:
%           min   |y - x*b|
%           s.t.  -M*z <=    b    <= M*z
%                 -M*z <= t - x*b <= M*z
%                 sum(z) < p
%       where
%           p: number of predictors (p <= m)
%           z_j: indicator if predictor j is active
%
% solution = linearRegression(x,y,p)
% solution = linearRegression(x,y,p,M,interceptFlag,predFlag,normFlag,opt_params)
%
%REQUIRED INPUTS
% x: Explanatory variables [n x m matrix]
% y: Response variables [n x 1 or n x m matrix]
%   n: number of samples
%   m: number of attributes
% p: Number of predictors
%
%OPTIONAL INPUTS
% M: ... (default = 1E3)
% interceptFlag: Include ('true') or exclude ('false) the intercept
%  (default = 'true')
% predFlag: Flag to use predictors to use other predictors ('true') or not
%  ('false'). If false, predictors can only predict responses.
%  (default = 'false')
% normFlag: Flag to normalize ('norm') or standardize ('stand') the data
%  (default = 'none')
%   norm: x_norm = (x - min(x))./(max(x) - min(x));
%   stand: x_stand = (x - mean(x))./(std(x));
% opt_params: Structure containing Gurobi parameters. A full list may be
%   found on the Parameter page of the Gurobi reference manual:
%       https://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters
%
%OUTPUT
% The solution structure will contains the following fields:
%   solution.B: Regression coefficients [(m+1) x 1 vector or (m+1) x m matrix]
%   solution.w: Loss [n x m matrix]
%   solution.z: Indicator if predictor j is active [m x 1 vector]
%               1: predictor, 0: not a predictor
%   solution.y: Estimated y values (X*B) [n x 1 or n x m matrix]
%   solution.status: Status (optimal, infeasible)

%% Check Inputs

if (nargin < 3)
    error('myfuns:linearRegression:NotEnoughInputs', ...
        'Not enough inputs: need x, y, and p');
else
    if ~ismatrix(x)
        error('myfuns:linearRegression:IncorrectType', ...
            '"x" needs to be a matrix');
    end
    if size(x,1) ~= size(y,1)
        error('myfuns:linearRegression:IncorrectSize', ...
            'x and y must have the same number of samples, n');
    end
end
[n,m] = size(x); % n observations, m features

if ~exist('M','var') || isempty(M)
    M = 1E3; % ...
end

if ~exist('interceptFlag','var') || isempty(interceptFlag)
    interceptFlag = 'true'; % include intercept
end

if ~exist('predFlag','var') || isempty(predFlag)
    predFlag = 'false'; % predictors can only predict responses
end

if ~exist('normFlag','var') || isempty(normFlag)
    normFlag = 'none'; % don't normalize or standardize the data
elseif strcmp(normFlag,'norm')
    max_x = nanmax(x); min_x = nanmin(x);
    x = (x - min_x)./(max_x - min_x);
    max_y = nanmax(y); min_y = nanmin(y);
    y = (y - min_y)./(max_y - min_y);
elseif strcmp(normFlag,'stand')
    avg_x = nanmean(x); stdev_x = nanstd(x);
    x = (x - avg_x)./stdev_x;
    avg_y = nanmean(y); stdev_y = nanstd(y);
    y = (y - avg_y)./stdev_y;
else
    error('myfuns:linearRegression:IncorrectInput', ...
            'normFlag must be "none", "norm", or "stand"');
end

if ~exist('opt_params','var')
    opt_params.FeasibilityTol = 1e-9; % all values must be satisfied to this tolerance (min value)
    opt_params.OutputFlag = 0; % silence gurobi
    opt_params.DisplayInterval = 1; % frequency at which log lines are printed (in seconds)
else
    if ~isfield(opt_params,'FeasibilityTol')
        opt_params.FeasibilityTol = 1e-9; % all values must be satisfied to this tolerance (min value)
    end
    if ~isfield(opt_params,'OutputFlag')
        opt_params.OutputFlag = 0; % silence gurobi
    end
    if ~isfield(opt_params,'DisplayInterval')
        opt_params.DisplayInterval = 1; % frequency at which log lines are printed (in seconds)
    end
end

%% Build Model

if strcmp(interceptFlag,'true')
    % add column of 1s for B0
    X = [ones(n,1), x];
else
    X = x;
end

if size(y,2) > 1 % multivariate linear regression
    model = build_multivariate(X,y,m,n,p,M,interceptFlag,predFlag);
else % simple linear regression
    model = build_simple(X,y,m,n,p,M,interceptFlag,predFlag);
end

%% Run Model

sol = gurobi(model,opt_params);

if size(y,2) > 1 % multivariate linear regression
    if strcmp(sol.status,'OPTIMAL')
        solution.objectiveValue = sol.objval;
        if strcmp(normFlag,'none')
            solution.B = reshape(sol.x(1:m*(m+1)),m+1,m);
        elseif strcmp(normFlag,'norm')
            x = x.*(max_x - min_x) + min_x;
            y = y.*(max_y - min_y) + min_y;
            X = [ones(n,1), x];
            temp_B = reshape(sol.x(1:m*(m+1)),m+1,m);
            solution.B = NaN(m+1,m);
            solution.B(1,:) = (max_y - min_y) .* (temp_B(1,:) - (min_x./(max_x - min_x))*temp_B(2:end,:)) + min_y;
            solution.B(2:end,:) = (cell2mat(arrayfun(@(x) bsxfun(@rdivide, max_y-min_y,max_x(x)-min_x(x)), 1:m, 'Uni',false)')) .* temp_B(2:end,:);
        elseif strcmp(normFlag,'stand')
            x = x.*stdev_x + avg_x;
            y = y.*stdev_y + avg_y;
            X = [ones(n,1), x];
            temp_B = reshape(sol.x(1:m*(m+1)),m+1,m);
            solution.B = NaN(m+1,m);
            solution.B(1,:) = stdev_y .* (temp_B(1,:) - (avg_x./stdev_x)*temp_B(2:end,:)) + avg_y;
            solution.B(2:end,:) = (cell2mat(arrayfun(@(x) bsxfun(@rdivide, stdev_y,stdev_x(x)), 1:m, 'Uni',false)')) .* temp_B(2:end,:);
        end
        solution.w = reshape(sol.x(m*(m+1)+(1:m*n)),n,m);
        solution.t = reshape(sol.x(m*(m+1)+m*n+(1:m*n)),n,m);
        solution.z = (sol.x((m*(m+1)+2*m*n)+(1:m)))';
        solution.y = X*solution.B;
        % 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
        solution.Rsquared = arrayfun(@(x) 1 - (sum((y(:,x) - solution.y(:,x)).^2) / sum((y(:,x) - mean(y(:,x))).^2)),1:m);
        % 1 - (sum((y - yhat)^2)/(n-p-1)) / (sum((y - mean(y))^2)/(n-1))
        solution.Rsquared_adjusted = arrayfun(@(x) 1 - (sum((y(:,x) - solution.y(:,x)).^2)/(n-p-1)) / (sum((y(:,x) - mean(y(:,x))).^2)/(n-1)),1:m);
        % (sum((yhat - mean(y))^2)/p) / (sum((y - yhat)^2)/(n-p-1))
        solution.overall_F = arrayfun(@(x) (sum((solution.y(:,x) - mean(y(:,x))).^2)/p) / (sum((y(:,x) - solution.y(:,x)).^2)/(n-p-1)),1:m);
        solution.overall_pvalue = fpdf(solution.overall_F,p,n-p-1);
        % ((sum((y - yred)^2) - sum((y - yhat)^2))/1) / (sum((y - yhat)^2)/(n-p-1))
        solution.F = NaN(m+1,m);
        for ii = 1:(m+1)
            Bred = solution.B; Bred(ii,:) = 0;
            yred = X*Bred;
            solution.F(ii,:) = arrayfun(@(x) ((sum((y(:,x) - yred(:,x)).^2) - sum((y(:,x) - solution.y(:,x)).^2))/1) / (sum((y(:,x) - solution.y(:,x)).^2)/(n-p-1)),1:m);
        end
        solution.pvalue = fpdf(solution.F,1,n-p);
    else
        solution.objectiveValue = NaN;
        solution.B = NaN(m+1,m);
        solution.w = NaN(n,m);
        solution.z = NaN(m,1);
        solution.y = NaN(n,m);
        solution.Rsquared = NaN(1,m);
        solution.Rsquared_adjusted = NaN(1,m);
        solution.F = NaN(1,m);
        solution.pvalue = NaN(1,m);
    end
else % simple linear regression
    if strcmp(sol.status,'OPTIMAL') && strcmp(interceptFlag,'true')
        solution.objectiveValue = sol.objval;
        if strcmp(normFlag,'none')
            solution.B = sol.x(1:(m+1));
        elseif strcmp(normFlag,'norm')
            x = x.*(max_x - min_x) + min_x;
            y = y.*(max_y - min_y) + min_y;
            X = [ones(n,1), x];
            temp_B = sol.x(1:(m+1));
            solution.B = NaN(m+1,1);
            solution.B(1,:) = (max_y - min_y) * (temp_B(1) - (min_x./(max_x - min_x))*temp_B(2:end)) + min_y;
            solution.B(2:end,:) = ((max_y-min_y)./(max_x-min_x)) .* temp_B(2:end)';
        elseif strcmp(normFlag,'stand')
            x = x.*stdev_x + avg_x;
            X = [ones(n,1), x];
            temp_B = sol.x(1:(m+1));
            solution.B = NaN(m+1,1);
            solution.B(1,:) = stdev_y*(temp_B(1) - (avg_x./stdev_x)*temp_B(2:end)) + avg_y;
            solution.B(2:end,:) = (stdev_y./stdev_x) .* temp_B(2:end)';
        end
        solution.w = sol.x((m+1)+(1:n));
        solution.z = sol.x((m+1+n)+(1:m));
        solution.y = X*solution.B;
        % 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
        solution.Rsquared = 1 - (sum((y - solution.y).^2) / sum((y - mean(y)).^2));
        % 1 - (sum((y - yhat)^2)/(n-p-1)) / (sum((y - mean(y))^2)/(n-1))
        solution.Rsquared_adjusted = 1 - ((sum((y - solution.y).^2)/(n-p-1)) / (sum((y - mean(y)).^2))/(n-1));
        % (sum((yhat - mean(y))^2)/p) / (sum((y - yhat)^2)/(n-p-1))
        solution.overall_F = ((sum((solution.y - mean(y)).^2)/(p)) / (sum((y -solution.y).^2))/(n-p-1));
        solution.overall_pvalue = fpdf(solution.overall_F,p,n-p-1);
        % ((sum((y - yred)^2) - sum((y - yhat)^2))/1) / (sum((y - yhat)^2)/(n-p-1))
        solution.F = NaN(m+1,1);
        for ii = 1:(m+1)
            Bred = solution.B; Bred(ii) = 0;
            yred = X*Bred;
            solution.F(ii) = ((sum((y - yred).^2) - sum((y - solution.y).^2))/1) / (sum((y - solution.y).^2)/(n-p-1));
        end
        solution.pvalue = fpdf(solution.F,1,n-p);
    elseif strcmp(sol.status,'OPTIMAL') && strcmp(interceptFlag,'false')
        solution.objectiveValue = sol.objval;
        if strcmp(normFlag,'none')
            solution.B = sol.x(1:m);
        elseif strcmp(normFlag,'norm')
            x = x.*(max_x - min_x) + min_x;
            y = y.*(max_y - min_y) + min_y;
            X = [ones(n,1), x];
            temp_B = sol.x(1:m);
            solution.B = ((max_y-min_y)./(max_x-min_x)) .* temp_B(2:end)';
        elseif strcmp(normFlag,'stand')
            x = x.*stdev_x + avg_x;
            X = [ones(n,1), x];
            temp_B = sol.x(1:(m+1));
            solution.B = (stdev_y./stdev_x) .* temp_B(2:end)';
        end
        solution.w = sol.x(m+(1:n));
        solution.z = sol.x((m+n)+(1:m));
        solution.y = X*solution.B;
        % 1 - sum((y - yhat)^2)/sum((y - mean(y))^2)
        solution.Rsquared = 1 - (sum((y - solution.y).^2) / sum((y - mean(y)).^2));
        % 1 - (sum((y - yhat)^2)/(n-p-1)) / (sum((y - mean(y))^2)/(n-1))
        solution.Rsquared_adjusted = 1 - ((sum((y - solution.y).^2)/(n-p-1)) / (sum((y - mean(y)).^2))/(n-1));
        % (sum((yhat - mean(y))^2)/p) / (sum((y - yhat)^2)/(n-p-1))
        solution.overall_F = ((sum((solution.y - mean(y)).^2)/(p)) / (sum((y -solution.y).^2))/(n-p-1));
        solution.overall_pvalue = fpdf(solution.overall_F,p,n-p-1);
        % ((sum((y - yred)^2) - sum((y - yhat)^2))/1) / (sum((y - yhat)^2)/(n-p-1))
        solution.F = NaN(m,1);
        for ii = 1:m
            Bred = solution.B; Bred(ii) = 0;
            yred = X*Bred;
            solution.F(ii) = ((sum((y - yred).^2) - sum((y - solution.y).^2))/1) / (sum((y - solution.y).^2)/(n-p-1));
        end
        solution.pvalue = fpdf(solution.F,1,n-p);
    else
        solution.objectiveValue = NaN;
        solution.B = NaN(m+1,1);
        solution.w = NaN(n,m);
        solution.z = NaN(m,1);
        solution.y = NaN(n,1);
        solution.Rsquared = NaN;
        solution.Rsquared_adjusted = NaN;
        solution.overall_F = NaN;
        solution.overall_pvalue = NaN;
        solution.F = NaN(m+1,1);
        solution.pvalue = NaN(m+1,1);
    end
end
solution.status = sol.status;

end

%% Simple Linear Regression

function model = build_simple(X,y,m,n,p,M,interceptFlag,predFlag)

if strcmp(interceptFlag,'true')
    B_Mz = [zeros(m,1), eye(m)];
    
    % A matrix
    A = [
        % |y - XB|
                  -X,    -eye(n), zeros(n,m); %  y - XB < w
                   X,    -eye(n), zeros(n,m); % -y + XB < w
        % -Mz < B < Mz
               -B_Mz, zeros(m,n), -M.*eye(m); % -Mz < B
                B_Mz, zeros(m,n), -M.*eye(m); % B < Mz
        % sum(z) < p
        zeros(1,m+1), zeros(1,n),  ones(1,m); % sum(z)<p
        ];
    
    % Right Hand Side
    rhs = [
        % |y - XB|
        -y; %  y - XB < w
         y; % -y + XB < w
        % -Mz < B < Mz
        zeros(2*m,1);
        % sum(z) < p
        p;
        ];

    % Lower Bound
    lb = [
        -M*ones(m+1,1); % B
            zeros(n,1); % w
            zeros(m,1); % z
        ];

    % Upper Bound
    ub = [
        M*ones(m+1+n,1); % B, w
              ones(m,1); % z
        ];

    % Objective
    obj = [
        zeros(m+1,1); % B
        ones(n,1)./n; % w
          zeros(m,1); % z
        ];

    % Sense on RHS
    sense = repmat('<',2*n+2*m+1,1);

    % Variable Types
    vtype = [
        repmat('C',m+1+n,1); % B, w
            repmat('B',m,1); % z
        ];
else
    B_Mz = eye(m);
    
    % A matrix
    A = [
        % |y - XB|
                  -X,    -eye(n), zeros(n,m); %  y - XB < w
                   X,    -eye(n), zeros(n,m); % -y + XB < w
        % -Mz < B < Mz
               -B_Mz, zeros(m,n), -M.*eye(m); % -Mz < B
                B_Mz, zeros(m,n), -M.*eye(m); % B < Mz
        % sum(z) < p
          zeros(1,m), zeros(1,n),  ones(1,m); % sum(z)<p
        ];
    
    % Right Hand Side
    rhs = [
        % |y - XB|
        -y; %  y - XB < w
         y; % -y + XB < w
        % -Mz < B < Mz
        zeros(2*m,1);
        % sum(z) < p
        p;
        ];

    % Lower Bound
    lb = [
          -M*ones(m,1); % B
            zeros(n,1); % w
            zeros(m,1); % z
        ];

    % Upper Bound
    ub = [
          M*ones(m+n,1); % B, w
              ones(m,1); % z
        ];

    % Objective
    obj = [
          zeros(m,1); % B
        ones(n,1)./n; % w
          zeros(m,1); % z
        ];

    % Sense on RHS
    sense = repmat('<',2*n+2*m+1,1);

    % Variable Types
    vtype = [
          repmat('C',m+n,1); % B, w
            repmat('B',m,1); % z
        ];
end

% Structure
model.A = sparse(A);
model.rhs = rhs;
model.lb = lb;
model.ub = ub;
model.obj = obj;
model.modelsense = 'min';
model.sense = sense;
model.vtype = vtype;

end

%% Multivariate Linear Regression

function model = build_multivariate(X,y,m,n,p,M,interceptFlag,predFlag)

% restructure X for A matrix
A_X = repmat({sparse(X)},1,m);
A_X = blkdiag(A_X{:});

if strcmp(interceptFlag,'true')
    % -Mz < B < Mz
    % +/- B
    A_B_r = repmat({sparse([zeros(m,1), eye(m)])},1,m);
    A_B_r = blkdiag(A_B_r{:});
    % - Mz
    A_Mz_r = repmat(-M.*eye(m),m,1);

    % -M(1-z) < B < M(1-z)
    % Mz
    A_Mz_p = repmat({M.*ones(m+1,1)},1,m);
    A_Mz_p = blkdiag(A_Mz_p{:});
else
    % -Mz < B < Mz
    % +/- B
    A_B_r = repmat({sparse(eye(m))},1,m);
    A_B_r = blkdiag(A_B_r{:});
    % - Mz
    A_Mz_r = repmat(-M.*eye(m),m,1);

    % -M(1-z) < B < M(1-z)
    % Mz
    A_Mz_p = repmat({M.*ones(m,1)},1,m);
    A_Mz_p = blkdiag(A_Mz_p{:});
end

% -Mz < t - x < Mz
A_tx = repmat({M.*ones(n,1)},m,1);
A_tx = blkdiag(A_tx{:});

if strcmp(predFlag,'true')
    % A matrix
    A = sparse([
        %                  B,                  w,                  t,              z;
        % |t - XB|
                        -A_X,          -eye(m*n),           eye(m*n),   zeros(m*n,m); %  t - XB < w
                         A_X,          -eye(m*n),          -eye(m*n),   zeros(m*n,m); % -t + XB < w
        % -Mz_k < B_kj < Mz_k
                      -A_B_r,     zeros(m*m,m*n),     zeros(m*m,m*n),         A_Mz_r;  % -Mz < B
                       A_B_r,     zeros(m*m,m*n),     zeros(m*m,m*n),         A_Mz_r;  % B < Mz
        % -Mz_k < t_ik - x_ik < Mz_k
          zeros(m*n,m*(m+1)),     zeros(m*n,m*n),           eye(m*n),          -A_tx; %  t - x < Mz
          zeros(m*n,m*(m+1)),     zeros(m*n,m*n),          -eye(m*n),          -A_tx; % -Mz < t - x
        % sum(z_k) < p
            zeros(1,m*(m+1)),       zeros(1,m*n),       zeros(1,m*n),      ones(1,m); % sum(z) < p
        ]);

    % Right Hand Side
    rhs = [
        % |t - XB|
        zeros(2*m*n,1);
        % -Mz < B < Mz
        zeros(2*m*m,1);
        % -Mz < t - x < Mz
         y(:);
        -y(:);
        % sum(z) < p
        p;
        ];    
else
    % A matrix
    A = sparse([
        %                  B,                  w,                  t,              z;
        % |t - XB|
                        -A_X,          -eye(m*n),           eye(m*n),   zeros(m*n,m); %  t - XB < w
                         A_X,          -eye(m*n),          -eye(m*n),   zeros(m*n,m); % -t + XB < w
        % -Mz_k < B_kj < Mz_k
                      -A_B_r,     zeros(m*m,m*n),     zeros(m*m,m*n),         A_Mz_r;  % -Mz < B
                       A_B_r,     zeros(m*m,m*n),     zeros(m*m,m*n),         A_Mz_r;  % B < Mz
        % - M(1-z_k) < B_jk < M(1-z_k)
               -eye(m*(m+1)), zeros(m*(m+1),m*n), zeros(m*(m+1),m*n),         A_Mz_p;  % -M(1-z) < B
                eye(m*(m+1)), zeros(m*(m+1),m*n), zeros(m*(m+1),m*n),         A_Mz_p;  % B < M(1-z)
        % -Mz_k < t_ik - x_ik < Mz_k
          zeros(m*n,m*(m+1)),     zeros(m*n,m*n),           eye(m*n),          -A_tx; %  t - x < Mz
          zeros(m*n,m*(m+1)),     zeros(m*n,m*n),          -eye(m*n),          -A_tx; % -Mz < t - x
        % sum(z_k) < p
            zeros(1,m*(m+1)),       zeros(1,m*n),       zeros(1,m*n),      ones(1,m); % sum(z) < p
        ]);

    % Right Hand Side
    rhs = [
        % |t - XB|
        zeros(2*m*n,1);
        % -Mz < B < Mz
        zeros(2*m*m,1);
        % -M(1-z) < B < M(1-z)
        M.*ones(2*m*(m+1),1);
        % -Mz < t - x < Mz
         y(:);
        -y(:);
        % sum(z) < p
        p;
        ];
end

% B lower/upper bounds
B_bounds = M*ones(m*(m+1),1); % B
B_bounds(2:(m+2):end) = 0;

% Lower Bound
lb = [
            -B_bounds; % B
         zeros(m*n,1); % w
%             -10.*y(:); % t
      -M.*ones(m*n,1); % t
           zeros(m,1); % z
    ];

% Upper Bound
ub = [
             B_bounds; % B
     Inf.*ones(m*n,1); % w
%              10.*y(:); % t
       M.*ones(m*n,1); % t
            ones(m,1); % z
    ];

% Objective
obj = [
    zeros(m*(m+1),1); % B
      ones(m*n,1)./n; % w
        zeros(m*n,1); % t
          zeros(m,1); % z
    ];

% Sense on RHS
sense = repmat('<',numel(rhs),1);

% Variable Types
vtype = [
     repmat('C',m*(m+1),1); % B
         repmat('C',m*n,1); % w
         repmat('C',m*n,1); % t
           repmat('B',m,1); % z
    ];

% Variable Names
[a,b] = ndgrid(0:m,1:m);
B = regexp(sprintf('B_%d_%d,',[a(:)'; b(:)']),',','split'); B(end) = [];
B = reshape(B,m+1,[]);
[a,b] = ndgrid(1:n,1:m);
w = regexp(sprintf('w_%d_%d,',[a(:)'; b(:)']),',','split'); w(end) = [];
w = reshape(w,n,[]);
t = regexp(sprintf('t_%d_%d,',[a(:)'; b(:)']),',','split'); t(end) = [];
t = reshape(t,n,[]);
z = regexp(sprintf('z_%d,',1:m),',','split'); z(end) = [];
z = z';

% Structure
model.A = A;
model.rhs = rhs;
model.lb = lb;
model.ub = ub;
model.obj = obj;
model.modelsense = 'min';
model.sense = sense;
model.vtype = vtype;
model.varname = [B(:); w(:); t(:); z];

end


