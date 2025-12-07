function PSTH = activation_identity(ConvTrace, Params)
% Params.k : scaling factor (spikes/sec per signal unit)
PSTH = 10000.* ConvTrace(1,:); % assume 1D
end