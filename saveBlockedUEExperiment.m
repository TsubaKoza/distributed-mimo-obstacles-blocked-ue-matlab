function saveBlockedUEExperiment(experiment)
%SAVEBLOCKEDUEEXPERIMENT Named MAT variables plus CSVs and a run summary.
c=experiment.config; folder=c.resultsDir;
if ~exist(folder,'dir'), mkdir(folder); end
save(fullfile(folder,'blocked_ue_experiment.mat'),'-struct','experiment','-v7.3');
fields={'channelRankingTable','singleAPMinimumPowerTable','allPowerCombinationTable', ...
    'feasiblePowerCombinationTable','topPowerCombinationTable','scenarioAttemptTable','targetCandidateTable'};
for i=1:numel(fields)
    writetable(experiment.(fields{i}),fullfile(folder,[fields{i},'.csv']));
end
if experiment.minimumPowerSolution.Found
    writetable(experiment.minimumPowerSolution.Combination,fullfile(folder,'minimumPowerCombination.csv'));
    writetable(experiment.minimumPowerSolution.APAllocation,fullfile(folder,'minimumPowerAPAllocation.csv'));
end
writetable(struct2table(experiment.correlationResults),fullfile(folder,'correlationResults.csv'));
if isfield(experiment,'continuousPowerReference') && experiment.continuousPowerReference.Feasible
    writetable(experiment.continuousPowerReference.APAllocation,fullfile(folder,'continuousPowerReference.csv'));
end
f=fopen(fullfile(folder,'RUN_SUMMARY.md'),'w','n','UTF-8'); assert(f>=0);
cleanup=onCleanup(@()fclose(f));
fprintf(f,'# Single blocked UE experiment: measured run\n\n');
fprintf(f,'MATLAB: %s. Accepted seed: %d; requested seed: %d.\n\n',version,c.seed,experiment.requestedConfig.seed);
fprintf(f,'APs: %d; UEs: %d; buildings: %d; serviceable blocked UEs: %d; target UE: %d; usable NLoS APs: %d.\n\n', ...
    c.M,c.K,numel(experiment.scenario.obstacles),nnz(experiment.classification.serviceableMask), ...
    experiment.targetBlockedUE.UEIndex,experiment.targetBlockedUE.NumNLoSAPs);
fprintf(f,'Selection: %s\n\n',experiment.targetBlockedUE.SelectionReason);
fprintf(f,'Carrier: %.6g GHz; bandwidth: %.6g MHz; target: %.6g Mbps; provisional NF: %.6g dB.\n\n', ...
    c.fc/1e9,c.bandwidth/1e6,c.targetRate/1e6,c.noiseFigure_dB);
fprintf(f,'Noise density: %.6g dBm/Hz; final noise: %.12g W (%.6g dBm). SINR = SNR (one served UE).\n\n', ...
    c.thermalNoiseDensity,experiment.noise.powerW,experiment.noise.powerDBm);
fprintf(f,'Reflection mode: %s; fixed magnitude: %g; phase: %g deg.\n\n', ...
    c.reflectionCoefficientMode,c.fixedReflectionMagnitude,c.fixedReflectionPhaseDeg);
fprintf(f,'Power step: %g W; cap: %g W/AP. Feasible combinations: %d / %d.\n\n', ...
    c.powerStep,c.maxTransmitPowerPerAP,height(experiment.feasiblePowerCombinationTable),height(experiment.allPowerCombinationTable));
if experiment.minimumPowerSolution.Found
    t=experiment.minimumPowerSolution.Combination;
    fprintf(f,'Minimum grid total power: %.12g W; SINR: %.9g (%.9g dB); rate: %.9g Mbps.\n\n', ...
        t.TotalPower_W,t.SINR_linear,t.SINR_dB,t.Rate_Mbps);
    fprintf(f,'|Rank|AP|Gain dB|Allocated W|Fraction|\n|---:|---:|---:|---:|---:|\n');
    a=experiment.minimumPowerSolution.APAllocation;
    for i=1:height(a)
        fprintf(f,'|%d|%d|%.9g|%.9g|%.9g|\n',a.Rank(i),a.APIndex(i),a.ChannelGain_dB(i),a.AllocatedPower_W(i),a.PowerFraction(i));
    end
end
fprintf(f,'\n## All usable NLoS links\n\n|Rank|AP|Gain dB|Reflection paths|Diffraction paths|Single AP required W (10 Mbps)|\n|---:|---:|---:|---:|---:|---:|\n');
t=experiment.channelRankingTable; t=t(t.FeasibleNLoSLink,:);
single=experiment.singleAPMinimumPowerTable;
for i=1:height(t)
    power=single.MinimumPowerFor10Mbps_W(single.APIndex==t.APIndex(i));
    fprintf(f,'|%d|%d|%.9g|%d|%d|%.9g|\n',t.Rank(i),t.APIndex(i),t.ChannelGain_dB(i),t.NumReflectionPaths(i),t.NumDiffractionPaths(i),power);
end
if isfield(experiment,'continuousPowerReference')
    ref=experiment.continuousPowerReference;
    fprintf(f,'\nContinuous reference (same APs): %.12g W; %.9g Mbps. This is separate from the discrete grid optimum.\n\n', ...
        ref.TotalPower_W,ref.Rate_Mbps);
end
r=experiment.correlationResults;
fprintf(f,'\nCorrelation: N=%d; Pearson=%.12g (p=%.12g); Spearman=%.12g (p=%.12g).\n\n%s\n\n%s\n\n', ...
    r.N,r.Pearson,r.PearsonPValue,r.Spearman,r.SpearmanPValue,r.PValueMethod,r.Interpretation);
fprintf(f,'Validation passed: %d; repeat geometry/channels/power table: %d.\n\n', ...
    experiment.validation.passed,experiment.validation.reproducibilityChecked);
if isfield(experiment.validation,'random3DRegressionChecked')
    fprintf(f,'Legacy random3D exact regression checked: %d; original Ray Tracing sanity: %d checks passed.\n\n', ...
        experiment.validation.random3DRegressionChecked,experiment.validation.legacySanity.numChecks);
end
fprintf(f,'Optimum scope: %s. See README_BLOCKED_UE.md for physical limitations and formula conventions.\n', ...
    experiment.minimumPowerSolution.Scope);
fprintf('Saved MAT, CSV, figures and RUN_SUMMARY.md in %s\n',folder);
end
