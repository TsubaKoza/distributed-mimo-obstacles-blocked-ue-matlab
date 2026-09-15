function reference = computeContinuousPowerReference(ranking,apIndices,cfg)
%COMPUTECONTINUOUSPOWERREFERENCE Analytic/KKT reference for GRID comparison.
% u=sqrt(p): minimize sum(u.^2), subject to a.'*u>=sqrt(gamma*noise),
% 0<=u<=sqrt(Pmax). KKT gives u_i=min(sqrt(Pmax),lambda*a_i).
% This is an additional reference, NEVER a replacement for grid enumeration.
noise=computeExperimentNoise(cfg);
[found,rows]=ismember(apIndices,ranking.APIndex); assert(all(found));
t=ranking(rows,:); a=sqrt(t.EffectiveGainLinear);
assert(all(t.FeasibleNLoSLink) && all(a>0));
threshold=sqrt(noise.requiredSINR*noise.powerW);
cap=sqrt(cfg.maxTransmitPowerPerAP);
reference=struct('Feasible',sum(cap*a)>=threshold,'Scope', ...
    "Continuous-power optimum for the same selected APs under ideal coherent MRT; separate from the grid solution", ...
    'APAllocation',table(),'TotalPower_W',Inf,'SINR_linear',0,'Rate_Mbps',0);
if ~reference.Feasible, return; end
% Normalize to avoid an unnecessarily large Lagrange multiplier.
b=a/max(a); target=min(threshold/max(a),sum(cap*b));
low=0; high=1;
for bracket=1:1024
    if sum(b.*min(cap,high*b))>=target, break; end
    high=2*high;
end
assert(isfinite(high) && sum(b.*min(cap,high*b))>=target,'Unable to bracket continuous reference.');
for iteration=1:100
    mid=(low+high)/2;
    if sum(b.*min(cap,mid*b))>=target, high=mid; else, low=mid; end
end
power=min(cap,high*b).^2;
reference.APAllocation=t(:,{'Rank','APIndex','ChannelGain_dB','EffectiveGain_dB'});
reference.APAllocation.AllocatedPower_W=power;
reference.APAllocation.PowerFraction=power/sum(power);
reference.TotalPower_W=sum(power);
reference.SINR_linear=(sum(a.*sqrt(power)))^2/noise.powerW;
reference.Rate_Mbps=cfg.bandwidth*log2(1+reference.SINR_linear)/1e6;
assert(all(power<=cfg.maxTransmitPowerPerAP*(1+1e-12)));
assert(abs(reference.Rate_Mbps-cfg.targetRate/1e6)<1e-8*max(1,cfg.targetRate/1e6));
end
