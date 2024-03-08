function [Results] = lmeEEG_TFCE(T_Obs, maxTFCE, e_loc, E_H)
% lmeEEG_TFCE: Function to permorm TFCE
%   The function (adapted from ept_TFCE.m) requires functions from https://github.com/Mensen/ept_TFCE-matlab
%   See https://doi.org/10.1016/j.neuroimage.2012.10.027

% [Inputs]
% - T_Obs: Map of actual t-values. The map should be a "Channels x Timepoints"
%   or "Channels x Frequencies x Timepoints" matrix
% - maxTFCE: Maximum absolute TFCE value from each permutated t-map
% - e_loc: Electrode locations file created using EEGLAB
%   (https://sccn.ucsd.edu/eeglab/index.php)
% - E_H: Parameters of E and H (e.g., [0.66 2])

% [Output]
% - Results Structure including:
%       -T_Obs
%       -TFCE enhanced T_Obs
%       -Maximum absolute TFCE value from each permutated t-map
%       -T_Obs p-values from the permutation distribution
%       -Significance T_Obs mask (alpha = .05)

%% -- Error Checking -- %%

% Check Location File for Number of Channels
if ~isequal(size(T_Obs,1), length(e_loc))
    error ('Number of channels in data does not equal that of locations file')
end



%% Calculate channel neighbours
ChN = ept_ChN2(e_loc);

%% TFCE transformation of observed t-values
if ndims(T_Obs) == 2
    TFCE_Obs = ept_mex_TFCE2D(T_Obs, ChN, E_H);
elseif ndims(T_Obs) == 3
    TFCE_Obs = ept_mex_TFCE3D(T_Obs, ChN, E_H);
end

%% Calculating the p value from the permutation distribution
Alpha = .05;
nPerm = length(maxTFCE);
maxTFCE = sort([maxTFCE;max(abs(TFCE_Obs(:)))]);
maxTFCEcrit = maxTFCE(round(nPerm*(1-Alpha)));
Mask = abs(TFCE_Obs)>=maxTFCEcrit;
P_Values = NaN(size(TFCE_Obs,1),size(TFCE_Obs,2));
for idx = 1:size(TFCE_Obs,1)
    for jdx = 1:size(TFCE_Obs,2)
        P_Values(idx,jdx) = sum(abs(TFCE_Obs(idx,jdx))<=maxTFCE)/(nPerm+1);
    end
end

%% Output structure
Results.Obs                 = T_Obs;
Results.TFCE_Obs            = TFCE_Obs;
Results.maxTFCE             = maxTFCE;
Results.P_Values            = P_Values;
Results.Mask                = Mask;

end