function paths = pruneRayPaths(paths,cfg)
%PRUNERAYPATHS Optional NLoS pruning. Geometric LoS is NEVER removed.
% Threshold is |complexGain|^2 in dB, before the array response.
if ~cfg.enablePathPruning || isempty(paths), return; end
validateattributes(cfg.minimumPathPowerDB,{'numeric'},{'scalar','real','nonnan'});
assert(cfg.maximumNLoSPathsPerLink>=0 && ...
    (isinf(cfg.maximumNLoSPathsPerLink) || mod(cfg.maximumNLoSPathsPerLink,1)==0));
isLoS = string({paths.type})=="LoS";
power = abs([paths.complexGain]).^2;
indices = find(~isLoS & 10*log10(power)>=cfg.minimumPathPowerDB);
[~,order] = sort(power(indices),'descend');
indices = indices(order(1:min(numel(order),cfg.maximumNLoSPathsPerLink)));
keep = isLoS; keep(indices) = true;
paths = paths(keep); paths = paths(:);
end
