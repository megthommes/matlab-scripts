function solution = linearRegression(x,y,M,opt_params)
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
% solution = linearRegression(x,y,M,opt_params)
%
%REQUIRED INPUTS
% x: n x m matrix of explanatory variables
% y: n x 1 or n x m matrix of response variables
%   n: number of samples
%   m: number of explanatory variables
%
%OPTIONAL INPUTS
% M: ... (default = 1E3)
% opt_params: Structure containing Gurobi parameters. A full list may be
%   found on the Parameter page of the Gurobi reference manual:
%       https://www.gurobi.com/documentation/8.1/refman/parameters.html#sec:Parameters
%
%OUTPUT
% The solution structure will contains the following fields:
%   solution.B: Regression coefficients [(m+1) x m matrix]
%   solution.w: Loss [n x m matrix]
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

[n,m] = size(x); % n observations, m explanatory variables

% add column of 1s for B0
X = [ones(n,1), x];

if size(y,2) > 1 % multivariate linear regression
    model = build_multivariate(X,y,m,n,M);
else % simple linear regression
    model = build_simple(X,y,m,n,M);
end

%% Run Model

sol = gurobi(model,opt_params);

if size(y,2) > 1 % multivariate linear regression
    if strcmp(sol.status,'OPTIMAL')
        solution.B = reshape(sol.x(1:m*(m+1)),m+1,m);
        solution.w = reshape(sol.x(m*(m+1)+(1:m*n)),n,m);
        solution.y = X*solution.B;
    else
        solution.B = NaN(m+1,m);
        solution.w = NaN(n,m);
        solution.y = NaN(n,m);
    end
else % simple linear regression
    if strcmp(sol.status,'OPTIMAL')
        solution.B = sol.x(1:(m+1));
        solution.w = sol.x((m+1)+(1:n));
        solution.y = X*solution.B;
    else
        solution.B = NaN(m+1,m);
        solution.w = NaN(n,m);
        solution.y = NaN(n,1);
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
       ones(n,1); % w
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
A_X = zeros(m*n,m*(m+1));
for ii = 1:m
    A_X( (1:n)+n*(ii-1) , (1:(m+1))+(m+1)*(ii-1) ) = X;
end

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
         ones(m*n,1); % w
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
