function solution = linearRegression(x,y,M,normFlag,opt_params)
%linearRegression Perform linear regression using optimization
%
% Perform linear regression of the form:
%           y = x*b + e
%       where
%           y: response/independent variables
%           x: explanatory/dependent variables
%           b: regression coefficients (to be solved for)
%           e: error/noise
% By solving a LP problem of the form:
%           min   |y - x*b|
%           s.t.  -M <= b <= M
%
% solution = linearRegression(x,y)
% solution = linearRegression(x,y,M,normFlag,opt_params)
%
%REQUIRED INPUTS
% x: Explanatory variables [n x m matrix]
% y: Response variables [n x 1 or n x m matrix]
%   n: number of samples
%   m: number of explanatory variables
%
%OPTIONAL INPUTS
% M: ... (default = 1E3)
% normFlag: Flag to normalize ('norm') or standardize ('stand') the data
%   (default = 'none')
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
%   solution.y: Estimated y values (X*B) [n x 1 or n x m matrix]
%   solution.status: Status (optimal, infeasible)
%   solution.Rsquared: (Multiple) Correlation coefficient [1 x 1 or 1 x m vector]

%% Check Inputs

if (nargin < 2)
    error('myfuns:linearRegression:NotEnoughInputs', ...
        'Not enough inputs: need x and y');
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

if ~exist('M','var') || isempty(M)
    M = 1E3; % ...
end

if ~exist('normFlag','var') || isempty(M)
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

%% Build Linear Programming Model

[n,m] = size(x); % n observations, m explanatory variables

% add column of 1s for B0
X = [ones(n,1), x];

if size(y,2) > 1 && ~isequal(x,y) % multivariate linear regression
    model = build_multivariate(X,y,m,n,M);
elseif size(y,2) > 1 && isequal(x,y) % multivariate linear regression
    model = build_multivariate_same(X,y,m,n,M);
else % simple linear regression
    model = build_simple(X,y,m,n,M);
end

%% Run Model

sol = gurobi(model,opt_params);

if size(y,2) > 1 % multivariate linear regression
    if strcmp(sol.status,'OPTIMAL')
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
        solution.y = X*solution.B;
        p = m;
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
        solution.B = NaN(m+1,m);
        solution.w = NaN(n,m);
        solution.y = NaN(n,m);
        solution.Rsquared = NaN(1,m);
        solution.Rsquared_adjusted = NaN(1,m);
        solution.F = NaN(m+1,m);
        solution.pvalue = NaN(m+1,m);
    end
else % simple linear regression
    if strcmp(sol.status,'OPTIMAL')
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
        solution.y = X*solution.B;
        p = m;
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
    else
        solution.B = NaN(m+1,1);
        solution.w = NaN(n,m);
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

function model = build_simple(X,y,m,n,M)

% A matrix
A = [
    % |y - XB|
    -X, -eye(n); %  y - XB < w
     X, -eye(n); % -y + XB < w
    ];

% Right Hand Side
rhs = [
    -y; %  y - XB < w
     y; % -y + XB < w
    ];

% Lower Bound
lb = [
    -M*ones(m+1,1); % B
        zeros(n,1); % w
    ];

% Upper Bound
ub =  M*ones(m+1+n,1); % B, w

% Objective
obj = [
    zeros(m+1,1); % B
    ones(n,1)./n; % w
    ];

% Sense on RHS
sense = repmat('<',2*n,1);

% Variable Types
vtype = repmat('C',m+1+n,1); % B, w

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

function model = build_multivariate(X,y,m,n,M)

% restructure X for A matrix
A_X = repmat({X},1,m);
A_X = blkdiag(A_X{:});

% A matrix
A = [
    % |y - XB|
    -A_X, -eye(m*n); %  y - XB < w
     A_X, -eye(m*n); % -y + XB < w
    ];

% Right Hand Side
rhs = [
    -y(:); %  y - XB < w
     y(:); % -y + XB < w
    ];

% Lower Bound
lb = [
    -M*ones(m*(m+1),1); % B
          zeros(m*n,1); % w
    ];

% Upper Bound
ub =  M*ones(m*(m+1)+m*n,1); % B, w

% Objective
obj = [
    zeros(m*(m+1),1); % B
      ones(m*n,1)./n; % w
    ];

% Sense on RHS
sense = repmat('<',2*m*n,1);

% Variable Types
vtype = repmat('C',m*(m+1)+m*n,1); % B, w

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

%% Multivariate Linear Regression (x==y)

function model = build_multivariate_same(X,y,m,n,M)

% restructure X for A matrix
A_X = repmat({X},1,m);
A_X = blkdiag(A_X{:});

% A matrix
A = [
    % |y - XB|
    -A_X, -eye(m*n); %  y - XB < w
     A_X, -eye(m*n); % -y + XB < w
    ];

% Right Hand Side
rhs = [
    -y(:); %  y - XB < w
     y(:); % -y + XB < w
    ];

% B lower/upper bounds
B_bounds = M*ones(m*(m+1),1); % B
B_bounds(2:(m+2):end) = 0;

% Lower Bound
lb = [
       -B_bounds; % B
    zeros(m*n,1); % w
    ];

% Upper Bound
ub = [
         B_bounds; % B
    M*ones(m*n,1); % w
    ];

% Objective
obj = [
    zeros(m*(m+1),1); % B
      ones(m*n,1)./n; % w
    ];

% Sense on RHS
sense = repmat('<',2*m*n,1);

% Variable Types
vtype = repmat('C',m*(m+1)+m*n,1); % B, w

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
