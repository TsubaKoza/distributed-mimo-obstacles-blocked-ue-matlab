function [allTable,feasibleTable,topTable,solution,singleTable] = sweepBlockedUEPower(ranking,amplitude,cfg)
%SWEEPBLOCKEDUEPOWER Exact enumeration on the chosen finite power grid.
noise=computeExperimentNoise(cfg); links=ranking(ranking.FeasibleNLoSLink,:);
singleTable=links(:,{'Rank','APIndex','ChannelGain_dB','EffectiveGain_dB'});
singleTable.MinimumPowerFor10Mbps_W=noise.requiredSINR*noise.powerW./links.EffectiveGainLinear;
singleTable.FeasibleWithin1W=singleTable.MinimumPowerFor10Mbps_W<=1;
singleTable.AchievedRateAt1W_Mbps=cfg.bandwidth*log2(1+links.EffectiveGainLinear/noise.powerW)/1e6;
% Explicit generic columns remain correct when users change QoS / AP cap.
singleTable.MinimumPowerForTargetRate_W=singleTable.MinimumPowerFor10Mbps_W;
singleTable.FeasibleWithinConfiguredCap=singleTable.MinimumPowerForTargetRate_W<=cfg.maxTransmitPowerPerAP;
singleTable.AchievedRateAtConfiguredCap_Mbps=cfg.bandwidth* ...
    log2(1+cfg.maxTransmitPowerPerAP*links.EffectiveGainLinear/noise.powerW)/1e6;
% Literal 10 Mbps columns must remain literal even with a changed targetRate.
gamma10=expm1(log(2)*10e6/cfg.bandwidth);
singleTable.MinimumPowerFor10Mbps_W=gamma10*noise.powerW./links.EffectiveGainLinear;
singleTable.FeasibleWithin1W=singleTable.MinimumPowerFor10Mbps_W<=1;
usable=links(links.EffectiveGainLinear>0,:);
selected=usable(1:min(height(usable),cfg.maxCandidateAPsForPowerSweep),:);
ids=selected.APIndex; n=numel(ids);
assert(n>0,'No usable NLoS APs.');
numLevels=ceil(cfg.maxTransmitPowerPerAP/cfg.powerStep)+1;
logCount=n*log(numLevels);
fprintf('Power sweep: %d APs, %d levels, %.0f combinations\n',n,numLevels,exp(logCount));
if logCount>log(cfg.maxPowerCombinations)+1e-12
    warning('BlockedUE:CombinationLimit','Requested sweep exceeds maxPowerCombinations=%g.',cfg.maxPowerCombinations);
    error('BlockedUE:SweepStopped','Reduce candidate APs or use a larger powerStep before allocating memory.');
end
levels=unique([0:cfg.powerStep:cfg.maxTransmitPowerPerAP,cfg.maxTransmitPowerPerAP]);
L=numel(levels); count=L^n; code=(0:count-1).'; power=zeros(count,n);
for j=1:n
    digits=mod(floor(code/L^(j-1)),L)+1;
    power(:,j)=reshape(levels(digits),[],1);
end
allTable=evaluateCoherentPowerAllocations(power,ids,amplitude(ids),cfg);
feasibleTable=sortrows(allTable(allTable.QoS_met,:), ...
    {'TotalPower_W','Rate_Mbps','CombinationID'},{'ascend','descend','ascend'});
topTable=feasibleTable(1:min(height(feasibleTable),cfg.numTopPowerSolutions),:);
solution=struct('Found',~isempty(feasibleTable),'CandidateAPIndices',ids, ...
    'Scope',"Discrete grid minimum over selected APs; not a continuous/global-all-AP optimum", ...
    'Combination',allTable([],:),'APAllocation',table());
if solution.Found
    solution.Combination=feasibleTable(1,:);
    allocated=power(feasibleTable.CombinationID(1),:).';
    allocation=selected(:,{'Rank','APIndex','ChannelGain_dB','EffectiveGain_dB'});
    allocation.AllocatedPower_W=allocated;
    allocation.PowerFraction=allocated/sum(allocated);
    solution.APAllocation=allocation;
end
end
