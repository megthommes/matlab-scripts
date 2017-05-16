function [FBA_solution] = FBA(model,exchRxns,minSumFlag)
%%FBA Solve a flux balance analysis problem
%
% Solves LP problems of the form: max/min c'*v
%                                 subject to S*v = b
%                                            lb <= v <= ub
% FBA_solution = FBA(model,exchRxns)
%
%REQUIRED INPUT
% The model structure must contain the following fields:
%   model.S:     Stoichiometric matrix
%   model.b:     Right hand side = dx/dt
%   model.c:     Objective coefficients
%   model.lb:    Lower bounds
%   model.ub:    Upper bounds
%
%OPTIONAL INPUTS
% exchRxns n x 3 matrix containing the lower (column 2) and upper (column 3)
%   bounds of exchange reactions at the specified indices (column 1)
% minSumFlag: true if would like to minimize the sum of fluxes (default=false)
%
%OUTPUT
% The FBA_solution structure will contains the following fields:
%   FBA_solution.objectiveValue    Objective value
%   FBA_solution.fluxes            Computed fluxes
%   FBA_solution.shadowPrices      Shadow prices (dual values) - for regular FBA only
%   FBA_solution.reducedCosts      Reduced costs - for regular FBA only
%   FBA_solution.status            Status (optimal, infeasible)
%
% Meghan Thommes 05/16/2017 - Added secondary opimization
% Meghan Thommes 4/9/2015 - Added vtype to FBA_model & renamed FBA_solution parameters
% Meghan Thommes 1/14/2015

%% Check Inputs

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

if exist('exchRxns','var') && ~isempty(exchRxns)
    if ~ismatrix(exchRxns)
        error('myfuns:FBA:IncorrectType', ...
            '"exchRxns" needs to be a matrix');
    elseif size(exchRxns,2) ~= 3
        error('myfuns:FBA:IncorrectType', ...
            '"exchRxns" needs to be an n x 3 matrix');
    else
        allExchRxns = full((sum(abs(model.S)==1,1) == 1) & (sum(model.S~=0) == 1))';
        allExchRxns = find(allExchRxns == 1);
        for i = 1:size(exchRxns,1)
            temp = find(allExchRxns == exchRxns(i,1),1);
            if isempty(temp)
                error('myfuns:FBA:IncorrectType', ...
                    '"exchRxns" contains a reaction that is not an exchange reaction');
            end
            clear temp
        end
    end
    clear allExchRxns
elseif isempty(exchRxns)
    clear exchRxns
end

if ~exist('minSumFlag','var') || isempty(minSumFlag)
    minSumFlag = false;
end
%%

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
FBA_params.DisplayInterval = 1; % frequency at which log lines are printed (in seconds)

% Uptake Exchange Reaction Bounds
if exist('exchRxns','var')
    FBA_model.lb(exchRxns(:,1)) = exchRxns(:,2);
    FBA_model.ub(exchRxns(:,1)) = exchRxns(:,3);
end

% Solving the Model: Linear Programming
solution = gurobi(FBA_model,FBA_params);

if strcmp(solution.status,'OPTIMAL')
    FBA_solution.objectiveValue = solution.objval;
    FBA_solution.fluxes = solution.x;
    FBA_solution.shadowPrices = solution.pi;
    FBA_solution.reducedCosts = solution.rc;
    
    % Second Linear Programming Problem (minimize total flux)
    if minSumFlag == 1
        % Build Irreversible Model
        pFBA_model.A = sparse([[model.S, -model.S]; [model.c; -model.c]']); % stoichiometric matrix, forward + reverse, with max growth objective function
        pFBA_model.obj = ones(size(pFBA_model.A,2),1); % linear objective vector, forward + reverse, all reactions
        pFBA_model.rhs = [model.b; FBA_model.obj'*solution.x]; % linear constraints, forward + reverse, with max growth objective function constraint
        % lower bounds on forward fluxes is max(0,lb), lower bound on reversed fluxes is max(0,-ub)
        pFBA_model.lb = [max(model.lb, zeros(numel(model.lb),1)); max(-model.ub, zeros(numel(model.ub),1))];
        % lower bounds on forward fluxes is max(0,ub), lower bound on reversed fluxes is max(0,-lb)
        pFBA_model.ub = [max(model.ub, zeros(numel(model.ub),1)); max(-model.lb, zeros(numel(model.lb),1))];
        pFBA_model.modelsense = 'min'; % minimize the objective function
        pFBA_model.sense = '='; % sense of the linear constraints for each row of A
        pFBA_model.vtype = 'C'; % continuous variables
        
        % Solve the Second Linear Programming Problem (minimize total flux)
        solution2 = gurobi(pFBA_model,FBA_params);
        
        % Minimzed Flux Values
        if strcmp(solution2.status,'OPTIMAL')
            FBA_solution.fluxes = solution2.x(1:end/2) - solution2.x(end/2+1:end); % combine forward and reverse fluxes
            
            % check objective value
            if (solution2.x(model.c==1) - solution.objval) > 1E-6
                FBA_solution.error = 'Minimization of absolute sum of fluxes failed';
            end
        else
            FBA_solution.error = 'Non-optimal solution when minimizing absolute sum of fluxes';
        end
    end
else
    FBA_solution.objectiveValue = NaN;
    
    % Second Linear Programming Problem (minimize total flux)
    if minSumFlag == 1
        FBA_solution.objectiveValue = NaN;
        FBA_solution.error = 'Non-optimal solution in initial optimization';
    end
end
FBA_solution.status = solution.status;




end