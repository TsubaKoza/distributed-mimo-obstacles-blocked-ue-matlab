function report = runBlockedUEPowerSanityChecks(baselineFile)
%RUNBLOCKEDUEPOWERSANITYCHECKS Analytic tests, legacy sanity, optional snapshot.
% Synthetic vectors below test algebra only; they never create research UEs.
legacy=runRayTracingSanityChecks();
c=configBlockedUEExperiment(struct('showFigures',false,'saveFigures',false,'saveResults',false));
noise=computeExperimentNoise(c);
assert(abs(noise.powerW-3.99052462993776e-13)<1e-25);
assert(abs(noise.requiredSINR-(sqrt(2)-1))<1e-14);
g=[1+2j,2-3j,0;-4+1j,1+4j,0]*1e-6;
[w,a,gain]=computeBlockedUEMRT(g);
assert(all(abs(gain-sum(abs(g).^2,1))<1e-25));
assert(all(w(:,3)==0) && a(3)==0);
assert(norm(w(:,1)-conj(g(:,1))/norm(g(:,1)))<1e-14);
assert(abs(g(:,1).'*w(:,1)-norm(g(:,1)))<1e-18);
p=[0 0;0.1 0;0 0.2;0.1 0.2];
t=evaluateCoherentPowerAllocations(p,[2;7],a(1:2),c);
assert(t.ReceivedDesiredPower_W(4)>sum(t.ReceivedDesiredPower_W(2:3)));
assert(~t.QoS_met(1) && t.Rate_Mbps(1)==0);
% Exact formula and required power are inverses (independent of sweep grid).
pmin=noise.requiredSINR*noise.powerW/gain(1);
exact=evaluateCoherentPowerAllocations(pmin,2,a(1),c);
assert(abs(exact.Rate_Mbps-c.targetRate/1e6)<1e-12);
% Pruning retains true geometric LoS even above all NLoS powers.
paths=struct('type',{"LoS","reflection","diffraction"}, ...
    'complexGain',{1e-9,1e-5,1e-8});
c.minimumPathPowerDB=-120;
kept=pruneRayPaths(paths,c);
assert(isequal(string({kept.type}),["LoS","reflection"]));
c.enablePathPruning=false; assert(isequal(pruneRayPaths(paths,c),paths));
% Single antenna case is a valid 1 x M matrix.
[w1,a1]=computeBlockedUEMRT([1j,2]);
assert(isequal(w1,[-1j,1]) && isequal(a1,[1,2]));
% Correlation handles no samples, one sample and a constant vector.
empty=table(zeros(0,1),zeros(0,1),'VariableNames',{'EffectiveGain_dB','MinimumPowerFor10Mbps_W'});
r=analyzeGainPowerCorrelation(empty); assert(r.N==0 && isnan(r.Pearson));
constant=table([1;1;1],[1;2;3],'VariableNames',empty.Properties.VariableNames);
r=analyzeGainPowerCorrelation(constant); assert(isnan(r.Pearson) && isnan(r.Spearman));
monotone=table([1;2;3],[3;2;1],'VariableNames',empty.Properties.VariableNames);
r=analyzeGainPowerCorrelation(monotone);
assert(abs(r.Spearman+1)<1e-14 && abs(r.SpearmanPValue-1/3)<1e-14);
% Fixed materials and non-overlapping outdoor geometry with nonzero yaw jitter.
c=configBlockedUEExperiment(struct('urbanYawJitterDeg',3));
rng(c.seed,'twister'); s=generateScenario3D(c); validateUrbanScenario3D(s,c);
% Continuous reference: analytic unconstrained case, active cap, infeasible.
toy=table([1;2],[1;2],[0;0],[0;0],[4;1],[true;true], ...
    'VariableNames',{'Rank','APIndex','ChannelGain_dB','EffectiveGain_dB','EffectiveGainLinear','FeasibleNLoSLink'});
c.bandwidth=1; c.thermalNoiseDensity=30; c.noiseFigure_dB=0; c.targetRate=1;
ref=computeContinuousPowerReference(toy,[1;2],c);
assert(abs(ref.TotalPower_W-0.2)<1e-12);
c.targetRate=log2(1+2.5^2);
ref=computeContinuousPowerReference(toy,[1;2],c);
assert(norm(ref.APAllocation.AllocatedPower_W-[1;0.25])<1e-12);
c.targetRate=log2(11); ref=computeContinuousPowerReference(toy,[1;2],c);
assert(~ref.Feasible);
regressionChecked=false;
if nargin>0 && ~isempty(baselineFile)
    saved=load(baselineFile,'baseline');
    c=configRayTracing3D('showFigures',false,'saveFigures',false,'saveResults',false);
    current=main_ray_tracing_3D(c);
    fields={'scenario','h_true','h_LoS','h_reflection','h_diffraction','pathInfo', ...
        'hasLoS','numReflections','numDiffractions','linkGainDB','outageMask'};
    for i=1:numel(fields)
        assert(isequaln(current.(fields{i}),saved.baseline.(fields{i})), ...
            'Legacy random3D differs from the pre-change snapshot: %s',fields{i});
    end
    regressionChecked=true;
end
report=struct('passed',true,'legacySanity',legacy,'random3DRegressionChecked',regressionChecked, ...
    'coverage',"Noise, MRT convention, coherent cross term, exact QoS power, pruning, zero gain, one antenna, correlation degeneracy/permutations, urban geometry with yaw, continuous optimum/cap/infeasibility, legacy sanity.");
fprintf('Blocked UE power sanity checks passed; legacy regression checked: %d\n',regressionChecked);
end
