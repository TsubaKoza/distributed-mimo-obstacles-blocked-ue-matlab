function tableOut = evaluateCoherentPowerAllocations(power,apIndices,amplitude,cfg)
%EVALUATECOHERENTPOWERALLOCATIONS One common symbol; amplitude sum BEFORE squaring.
noise=computeExperimentNoise(cfg);
received=abs(sqrt(power)*amplitude(:)).^2;
sinr=received/noise.powerW;
rate=cfg.bandwidth*log2(1+sinr)/1e6;
n=size(power,1);
tableOut=table((1:n).','VariableNames',{'CombinationID'});
for j=1:numel(apIndices)
    tableOut.(sprintf('P_AP%d_W',apIndices(j)))=power(:,j);
end
tableOut.TotalPower_W=sum(power,2);
tableOut.ReceivedDesiredPower_W=received;
tableOut.NoisePower_W=repmat(noise.powerW,n,1);
tableOut.SINR_linear=sinr;
tableOut.SINR_dB=10*log10(sinr);
tableOut.Rate_Mbps=rate;
tableOut.QoS_met=rate>=cfg.targetRate/1e6;
tableOut.NumActiveAPs=sum(power>0,2);
end
