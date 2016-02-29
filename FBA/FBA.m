function [FBA_solution] = FBA(model,excRxns)
%%FBA Solve a flux balance analysis problem
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% FBA_solution = FBA(model,excRxns)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.b:     Right hand side = dx/dt
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%   model.ub:    Upper bounds
%
%OPTIONAL INPUT
% excRxns n x 3 matrix containing the lower (column 2) and upper (column 3)
%   bounds of exchange reactions at the specified indices (column 1)
%
%OUTPUT
% The FBA_solution structure will always contain the following fields:
%   FBA_solution.objectiveValue    Objective value
%   FBA_solution.fluxes            Computed fluxes
%   FBA_solution.shadowPrices      Shadow prices (dual values)
%   FBA_solution.reducedCosts      Reduced costs
%   FBA_solution.status            Status (optimal, infeasible)
%
% Meghan Thommes 4/9/2015 - Added vtype to FBA_model & renamed FBA_solution parameters
% Meghan Thommes 1/14/2015

%%
% Check Inputs
if (nargin < 1)
    error('myfuns:FBA:NotEnoughInputs', ...
        'Not enough inputs: need a model');
else
    if ~isstruct(model)
        error('myfuns:FBA:IncorrectType', ...
            '"model" needs to be a structure');
    elseif ~isfield(model,'S') || ~isfield(model,'b') || ~isfield(model,'c') || ~isfield(model,'lb') || ~isfield(model,'ub')
        error('myfuns:FBA:IncorrectType', ...
            '"model" needs "S", "b", "c", "lb", and "ub" fields');
    end
end

if (nargin > 1)
    if ~ismatrix(excRxns)
        error('myfuns:FBA:IncorrectType', ...
            '"excRxns" needs to be a matrix');
    elseif size(excRxns,2) ~= 3
        error('myfuns:FBA:IncorrectType', ...
            '"excRxns" needs to be an n x 3 matrix');
    else
        allExcRxns = full((sum(abs(model.S)==1,1) == 1) & (sum(model.S~=0) == 1))';
        allExcRxns = find(allExcRxns == 1);
        for i = 1:size(excRxns,1)
            temp = find(allExcRxns == excRxns(i,1),1);
            if isempty(temp)
                error('myfuns:FBA:IncorrectType', ...
                    '"excRxns" contains a reaction that is not an exchange reaction');
            end
            clear temp
        end
    end
    clear allExcRxns
end

% Build Model
FBA_model.A = model.S; % stoichiometric matrix: linear constraint matrix [sparse matrix, rxns x mets]
FBA_model.obj = model.c; % linear objective vector for each each col of A (rxn in S) [dense vector]
FBA_model.rhs = model.b; % right-hand side vector for the linear constraints for each row of A (met in S) [dense vector]
FBA_model.lb = model.lb; % lower bounds for each col of A (rxn in S) [dense vector]
FBA_model.ub = model.ub; % upper bounds for each col of A (rxn in S) [dense vector]
FBA_model.modelsense = 'max'; % maximize objective function
FBA_model.sense = '='; % sense of the linear constraints for each row of A (met in S) [char array]
FBA_model.vtype = 'C'; % continuous variables

% Specify FBA Parameters
FBA_params.FeasibilityTol = 1e-9; % all values must be satisfied to this tolerance (min value)
FBA_params.OutputFlag = 0; % silence gurobi

% Uptake Exchange Reaction Bounds
if exist('excRxns','var')
    FBA_model.lb(excRxns(:,1)) = excRxns(:,2);
    FBA_model.ub(excRxns(:,1)) = excRxns(:,3);
end

% Solving the Model: Linear Programming
solution = gurobi(FBA_model,FBA_params);

if ~strcmp(solution.status,'INFEASIBLE')
    FBA_solution.objectiveValue = solution.objval;
    FBA_solution.fluxes = solution.x;
    FBA_solution.shadowPrices = solution.pi;
    FBA_solution.reducedCosts = solution.rc;
else
    FBA_solution.objectiveValue = NaN;
end
FBA_solution.status = solution.status;


end