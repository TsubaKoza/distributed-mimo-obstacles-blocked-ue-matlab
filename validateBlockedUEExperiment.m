function report = validateBlockedUEExperiment(e,reproduce)
%VALIDATEBLOCKEDUEEXPERIMENT Assertions on physics bookkeeping and experiment.
if nargin<2, reproduce=false; end
c=e.config; r=e.rayTracingResult; k=e.targetBlockedUE.UEIndex;
assert(all(~r.hasLoS(:,k)) && all(r.h_LoS(:,:,k)==0,'all'));
ids=e.minimumPowerSolution.CandidateAPIndices;
for m=ids(:).'
    types=string({r.pathInfo(m,k).paths.type});
    assert(any(types=="reflection" | types=="diffraction"));
    assert(~r.hasLoS(m,k) && ~any(types=="LoS") && ~r.outageMask(m,k));
end
t=e.allPowerCombinationTable;
p=t{:,2:1+numel(ids)};
assert(all(p>=0 & p<=c.maxTransmitPowerPerAP,'all'));
assert(max(abs(t.TotalPower_W-sum(p,2)))<1e-12);
assert(isequal(t.QoS_met,t.Rate_Mbps>=c.targetRate/1e6));
zero=all(p==0,2); assert(nnz(zero)==1 && ~any(t.QoS_met(zero)));
assert(t.ReceivedDesiredPower_W(zero)==0 && t.Rate_Mbps(zero)==0);
% Direct antenna-vector calculation, separate from the sweep's matrix sum.
rows=unique(round(linspace(1,height(t),min(19,height(t)))));
for row=rows
    signal=0;
    for j=1:numel(ids)
        m=ids(j); g=r.h_true(:,m,k); w=e.MRTWeights(:,m);
        hDL=conj(g);
        assert(abs(g.'*w-hDL'*w)<=1e-12*norm(g));
        signal=signal+sqrt(p(row,j))*(g.'*w);
    end
    expected=abs(signal)^2;
    assert(abs(expected-t.ReceivedDesiredPower_W(row))<=1e-11*max(expected,realmin));
end
for j=1:numel(ids)
    m=ids(j); single=e.singleAPMinimumPowerTable;
    required=single.MinimumPowerForTargetRate_W(single.APIndex==m);
    solo=all(p(:,setdiff(1:numel(ids),j))==0,2);
    achieved=p(solo & t.QoS_met,j);
    if required<=c.maxTransmitPowerPerAP*(1-1e-12)
        assert(~isempty(achieved));
        assert(min(achieved)>=required-1e-11);
        assert(min(achieved)<=required+c.powerStep+1e-11);
    elseif required>c.maxTransmitPowerPerAP*(1+1e-12)
        assert(isempty(achieved));
    end
end
assert(issorted(e.feasiblePowerCombinationTable.TotalPower_W));
if e.minimumPowerSolution.Found
    assert(e.minimumPowerSolution.Combination.TotalPower_W==min(t.TotalPower_W(t.QoS_met)));
    if isfield(e,'continuousPowerReference')
        assert(e.continuousPowerReference.TotalPower_W<=e.minimumPowerSolution.Combination.TotalPower_W+1e-10);
    end
end
if strcmpi(c.scenarioMode,"urbanStreetCanyon"), validateUrbanScenario3D(e.scenario,c); end
reproduced=false;
if reproduce
    fprintf('Reproducibility: re-generating geometry and ray tracing with accepted seed %d...\n',c.seed);
    state=rng; cleanup=onCleanup(@()rng(state));
    rng(c.seed,'twister'); s=generateScenario3D(c);
    assert(isequaln(s,e.scenario),'Same seed did not reproduce scenario.');
    [h,paths]=generateRayTracingChannels3D(s,c);
    assert(isequaln(h,r.h_true) && isequaln(paths,r.pathInfo), ...
        'Same seed did not reproduce ray-tracing channels/paths.');
    [~,a]=computeBlockedUEMRT(h(:,:,k));
    repeat=evaluateCoherentPowerAllocations(p,ids,a(ids),c);
    assert(isequaln(repeat,t),'Power table was not reproducible.');
    reproduced=true;
end
report=struct('passed',true,'reproducibilityChecked',reproduced, ...
    'message',"Blocked UE, retained paths, MRT conjugation, coherent power, cap, sums, QoS, zero-power, single-AP/grid agreement and geometry checks passed.");
fprintf('Experiment validation passed. Reproducibility checked: %d\n',reproduced);
end
