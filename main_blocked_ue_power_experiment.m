function experiment = main_blocked_ue_power_experiment(updates)
%MAIN_BLOCKED_UE_POWER_EXPERIMENT One blocked UE, coherent AP power sweep.
if nargin<1, updates=struct(); end
cfg=configBlockedUEExperiment(updates); requestedConfig=cfg;
noise=computeExperimentNoise(cfg);
fprintf('Reflection coefficient mode: %s\n',cfg.reflectionCoefficientMode);
if strcmpi(cfg.reflectionCoefficientMode,"fixed")
    fprintf('Reflection magnitude: %g\nReflection phase: %g deg\n', ...
        cfg.fixedReflectionMagnitude,cfg.fixedReflectionPhaseDeg);
end
fprintf('Carrier %.3g GHz; B %.3g MHz; NF %.3g dB (provisional); noise %.6g W (%.3f dBm)\n', ...
    cfg.fc/1e9,cfg.bandwidth/1e6,cfg.noiseFigure_dB,noise.powerW,noise.powerDBm);
fprintf('Target %.3g Mbps requires SINR %.6g (%.3f dB). Single UE: SINR = SNR.\n', ...
    cfg.targetRate/1e6,noise.requiredSINR,10*log10(noise.requiredSINR));
attempts=zeros(0,5); accepted=false; target=[];
for attempt=1:cfg.maxScenarioGenerationAttempts
    cfg.seed=requestedConfig.seed+attempt-1;
    rtCfg=cfg; rtCfg.showFigures=false; rtCfg.saveFigures=false; rtCfg.saveResults=false;
    fprintf('\nScenario attempt %d/%d, seed %d\n',attempt,cfg.maxScenarioGenerationAttempts,cfg.seed);
    result=main_ray_tracing_3D(rtCfg);
    classification=identifyServiceableBlockedUEs(result);
    [target,targetCandidateTable]=selectTargetBlockedUE(result,classification,cfg);
    blockedCount=nnz(classification.serviceableMask);
    attempts(end+1,:)=[attempt,cfg.seed,blockedCount, ...
        nnz(classification.completeOutageMask),nnz(targetCandidateTable.EligibleTarget)]; %#ok<AGROW>
    fprintf('Buildings=%d; all-AP blocked=%d; serviceable blocked=%d; eligible targets=%d\n', ...
        numel(result.obstacles),nnz(classification.allAPBlockedMask),blockedCount,nnz(targetCandidateTable.EligibleTarget));
    if blockedCount>=cfg.minBlockedUECount && ~isempty(target)
        accepted=true; break;
    end
end
attemptTable=array2table(attempts,'VariableNames', ...
    {'Attempt','Seed','ServiceableBlockedUECount','CompleteOutageUECount','EligibleTargetCount'});
if ~accepted
    if cfg.saveResults
        if ~exist(cfg.resultsDir,'dir'), mkdir(cfg.resultsDir); end
        writetable(attemptTable,fullfile(cfg.resultsDir,'scenario_attempts.csv'));
        save(fullfile(cfg.resultsDir,'scenario_generation_failure.mat'), ...
            'result','cfg','requestedConfig','classification','targetCandidateTable','attemptTable','-v7.3');
    end
    error('BlockedUE:NoTarget', ...
        'No accepted scenario after %d attempts. Last serviceable count=%d; requested=%d. No LoS/path was altered. Inspect attempt diagnostics.', ...
        cfg.maxScenarioGenerationAttempts,blockedCount,cfg.minBlockedUECount);
end
k=target.UEIndex;
fprintf('\nSelected blocked UE %d: %s\nNLoS AP count=%d; maximum sweep rate=%.6g Mbps\n', ...
    k,target.SelectionReason,target.NumNLoSAPs,target.MaximumSweepRate_Mbps);
[ranking,weights,amplitude]=buildBlockedUEChannelRanking(result,k);
[allTable,feasibleTable,topTable,solution,singleTable]=sweepBlockedUEPower(ranking,amplitude,cfg);
correlations=analyzeGainPowerCorrelation(singleTable);
continuousReference=computeContinuousPowerReference(ranking,solution.CandidateAPIndices,cfg);
experiment=struct('scenario',result.scenario,'targetBlockedUE',target, ...
    'channelRankingTable',ranking,'singleAPMinimumPowerTable',singleTable, ...
    'allPowerCombinationTable',allTable,'feasiblePowerCombinationTable',feasibleTable, ...
    'topPowerCombinationTable',topTable,'minimumPowerSolution',solution, ...
    'correlationResults',correlations,'config',cfg,'requestedConfig',requestedConfig, ...
    'rayTracingResult',result,'classification',classification,'MRTWeights',weights, ...
    'effectiveAmplitudes',amplitude,'noise',noise,'scenarioAttemptTable',attemptTable, ...
    'targetCandidateTable',targetCandidateTable,'continuousPowerReference',continuousReference);
experiment.validation=validateBlockedUEExperiment(experiment,cfg.verifyReproducibility);
fprintf('\n=== Channel Ranking ===\n'); disp(ranking);
fprintf('\n=== Single AP Minimum Power ===\n'); disp(singleTable);
fprintf('\nQoS feasible combinations: %d / %d\n',height(feasibleTable),height(allTable));
fprintf('\n=== Minimum Total Power Solution ===\n');
if solution.Found
    disp(solution.Combination); disp(solution.APAllocation);
else
    fprintf('No feasible grid combination.\n');
end
fprintf('\n=== Continuous Power Reference (same APs; separate from grid minimum) ===\n');
fprintf('Total %.9g W; rate %.9g Mbps\n',continuousReference.TotalPower_W,continuousReference.Rate_Mbps);
disp(continuousReference.APAllocation);
if cfg.showFigures || cfg.saveFigures
    plotBlockedUEExperiment(experiment);
end
if cfg.saveResults, saveBlockedUEExperiment(experiment); end
end
