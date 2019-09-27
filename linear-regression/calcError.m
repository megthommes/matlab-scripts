function [E] = calcError(y_meas,y_pred,errorFlag)
%calcError Calculates the mean error
%
% E = calcError(y_meas,y_pred,errorFlag)
%
%REQUIRED INPUTS
% y_meas: Measured values (n x m matrix)
% y_pred: Predicted values (n x m matrix)
% errorFlag: Type of error to calculate
%   MAE: mean absolute error
%   MSE: mean squared error
%   RMSE: root mean squared error
%
%OUTPUT
% E: Calculated error (1 x m vector)
%
% Meghan Thommes

%% Check Inputs

if (nargin < 3)
    error('myfuns:calcError:NotEnoughInputs', ...
        'Not enough inputs: need y_meas, y_pred, and errorFlag');
elseif ~strcmp(errorFlag,'MAE') && ~strcmp(errorFlag,'MSE') && ~strcmp(errorFlag,'RMSE')
    error('myfuns:calcError:IncorrectInput', ...
        'errorFlag must be "MAE", "MSE", or "RMSE"');
elseif size(y_meas,1) ~= size(y_pred,1) || size(y_meas,2) ~= size(y_pred,2)
    error('myfuns:calcError:IncorrectInput', ...
        'y_meas and y_pred must be the same size')
end

%% Calculate Error

[n,~] = size(y_meas);

if strcmp(errorFlag,'MAE') % Mean Absolute Error
    E = nansum(abs(y_meas - y_pred))./n;
elseif strcmp(errorFlag,'MSE') % Mean Squared Error
    E = nansum((y_meas - y_pred).^2)./n;
elseif strcmp(errorFlag,'RMSE') % Root Mean Squared Error
    E = sqrt(nansum((y_meas - y_pred).^2)./n);
end


end

