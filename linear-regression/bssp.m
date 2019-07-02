function solution = bssp(x,y,p,M,opt_params)
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
%           s.t.  -M*z <= b <= M*z
%                 sum(z) < p
%       where
%           p: number of predictors (p <= m)
%           z_j: indicator if predictor j is active
%
% solution = linearRegression(x,y)
% solution = linearRegression(x,y,p,M,opt_params)
%
%REQUIRED INPUTS
% x: Explanatory variables [n x m matrix]
% y: Response variables [n x 1 or n x m matrix]
%   n: number of samples
%   m: number of explanatory variables
%
%OPTIONAL INPUTS
% p: Number of predictors (default = m)
% M: ... (default = 1E3)
% opt_params: Structure containing Gurobi parameters. A full list may be
%   found on the Parameter page of the Gurobi reference manual:
%       https://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters
%
%OUTPUT
% The solution structure will contains the following fields:
%   solution.B: Regression coefficients [(m+1) x 1 vector or (m+1) x m matrix]
%   solution.w: Loss [n x m matrix]
%   solution.z: Indicator if predictor j is active [(m+1) x 1 vector]
%               1: predictor, 0: not a predictor (and corresponding B=0)
%   solution.y: Estimated y values (X*B) [n x 1 or n x m matrix]
%   solution.status: Status (optimal, infeasible)

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
[n,m] = size(x); % n observations, m explanatory variables

if ~exist('p','var') || isempty(p)
    p = m; % number of predictors
end

if ~exist('M','var') || isempty(M)
    M = 1E3; % ...
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

% add column of 1s for B0
X = [ones(n,1), x];

if size(y,2) > 1 % multivariate linear regression
    model = build_multivariate(X,y,m,n,p,M);
else % simple linear regression
    model = build_simple(X,y,m,n,p,M);
end

%% Run Model

sol = gurobi(model,opt_params);

if size(y,2) > 1 % multivariate linear regression
    if strcmp(sol.status,'OPTIMAL')
        solution.B = reshape(sol.x(1:m*(m+1)),m+1,m);
        solution.w = reshape(sol.x(m*(m+1)+(1:m*n)),n,m);
        solution.z = [1; sol.x((m*(m+1)+m*n)+(1:m))];
        solution.y = X*solution.B;
    else
        solution.B = NaN(m+1,m);
        solution.w = NaN(n,m);
        solution.z = NaN(m+1,1);
        solution.y = NaN(n,m);
    end
else % simple linear regression
    if strcmp(sol.status,'OPTIMAL')
        solution.B = sol.x(1:(m+1));
        solution.w = sol.x((m+1)+(1:n));
        solution.z = [1; sol.x((m+1+n)+(1:m))];
        solution.y = X*solution.B;
    else
        solution.B = NaN(m+1,1);
        solution.w = NaN(n,m);
        solution.z = NaN(m+1,1);
        solution.y = NaN(n,1);
    end
end
solution.status = sol.status;

end

%% Simple Linear Regression

function model = build_simple(X,y,m,n,p,M)

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
       ones(n,1); % w
      zeros(m,1); % z
    ];

% Sense on RHS
sense = repmat('<',2*n+2*m+1,1);

% Variable Types
vtype = [
    repmat('C',m+1+n,1); % B, w
        repmat('B',m,1); % z
    ];

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

function model = build_multivariate(X,y,m,n,p,M)

% restructure X for A matrix
A_X = repmat({X},1,m);
A_X = blkdiag(A_X{:});

% +/- B
A_B = repmat({[zeros(m,1), eye(m)]},1,m);
A_B = blkdiag(A_B{:});
% - Mz
A_Mz = repmat(-M.*eye(m),m,1);

% A matrix
A = [
    % |y - XB|
                    -A_X,      -eye(m*n),   zeros(m*n,m); %  y - XB < w
                     A_X,      -eye(m*n),   zeros(m*n,m); % -y + XB < w
    % -Mz < B < Mz
                    -A_B, zeros(m*m,m*n),           A_Mz;  % -Mz < B
                     A_B, zeros(m*m,m*n),           A_Mz;  % B < Mz
    % sum(z) < p
        zeros(1,m*(m+1)),   zeros(1,m*n),      ones(1,m); % sum(z) < p
    ];

% Right Hand Side
rhs = [
    % |y - XB|
    -y(:); %  y - XB < w
     y(:); % -y + XB < w
    % -Mz < B < Mz
    zeros(2*m*m,1);
    % sum(z) < p
    p;
    ];

% Lower Bound
lb = [
    -M*ones(m*(m+1),1); % B
          zeros(m*n,1); % w
            zeros(m,1); % z
    ];

% Upper Bound
ub = [
    M*ones(m*(m+1)+m*n,1); % B, w
                ones(m,1); % z
    ];

% Objective
obj = [
    zeros(m*(m+1),1); % B
         ones(m*n,1); % w
          zeros(m,1); % z
    ];

% Sense on RHS
sense = repmat('<',2*m*n+2*m*m+1,1);

% Variable Types
vtype = [
    repmat('C',m*(m+1)+m*n,1); % B, w
              repmat('B',m,1); % z
    ];

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
