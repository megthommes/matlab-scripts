function weight = getFluxWeight_mt(fluxValues,fluxScale)
%GETFLUXWEIGHT_MT Obtains color 
%
% weight = GETFLUXWEIGHT_MT(fluxValues,fluxScale)
%
%INPUT
% fluxValues:   Vector containing values
%
%OPTIONAL INPUT
% fluxScale:    Min/max flux value (default=[-1000,1000])
%
%OUTPUT
% weight        vector containing weight values for each input value
%
% Meghan Thommes 2/17/2017

% +/-Inf: 25
% +/-0.01: 8
% 0: 4

if nargin <2
    fluxScale = [-1000, 1000];
end

sample_flux = [0,0.01,fluxScale(2)]';
sample_weight = [4, 8, 25];
weight = interp1(sample_flux,sample_weight,abs(fluxValues));

end