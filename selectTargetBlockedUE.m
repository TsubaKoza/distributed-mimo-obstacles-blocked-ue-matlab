function [target,candidateTable] = selectTargetBlockedUE(result,classification,cfg)
%SELECTTARGETBLOCKEDUE Selection is based on post-RT eligibility and QoS.
noise=computeExperimentNoise(cfg);
ids=classification.serviceableIndices(:); count=zeros(numel(ids),1); maxRate=count;
for i=1:numel(ids)
    t=buildBlockedUEChannelRanking(result,ids(i));
    t=t(t.FeasibleNLoSLink & t.EffectiveGainLinear>0,:); count(i)=height(t);
    gain=t.EffectiveGainLinear(1:min(height(t),cfg.maxCandidateAPsForPowerSweep));
    power=(sum(sqrt(cfg.maxTransmitPowerPerAP*gain)))^2;
    maxRate(i)=cfg.bandwidth*log2(1+power/noise.powerW)/1e6;
end
eligible=count>=cfg.minimumTargetNLoSAPs & maxRate>=cfg.targetRate/1e6;
candidateTable=table(ids,count,maxRate,eligible,'VariableNames', ...
    {'UEIndex','UsableNLoSAPCount','MaximumSweepRate_Mbps','EligibleTarget'});
target=[];
if ~isempty(cfg.targetBlockedUEIndex)
    validateattributes(cfg.targetBlockedUEIndex,{'numeric'},{'scalar','integer','>=',1,'<=',cfg.K});
    row=find(ids==cfg.targetBlockedUEIndex & eligible,1);
    if isempty(row), return; end
    reason="Explicit UE index; all geometric, NLoS and maximum-power QoS checks passed.";
else
    candidates=find(eligible);
    if isempty(candidates), return; end
    [~,order]=sortrows([-count(candidates),-maxRate(candidates),ids(candidates)],[1 2 3]);
    row=candidates(order(1));
    reason="Most usable NLoS APs; ties resolved by maximum sweep rate, then UE index.";
end
target=struct('UEIndex',ids(row),'Position',result.UEpos(ids(row),:), ...
    'NumNLoSAPs',count(row),'MaximumSweepRate_Mbps',maxRate(row), ...
    'SelectionReason',reason);
end
