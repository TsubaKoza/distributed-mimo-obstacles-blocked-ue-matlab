function noise = computeExperimentNoise(cfg)
%COMPUTEEXPERIMENTNOISE Thermal density [dBm/Hz], B [Hz], receiver NF [dB].
dbm=cfg.thermalNoiseDensity+10*log10(cfg.bandwidth)+cfg.noiseFigure_dB;
noise=struct('powerDBm',dbm,'powerW',10^((dbm-30)/10), ...
    'requiredSINR',expm1(log(2)*cfg.targetRate/cfg.bandwidth));
assert(isfinite(noise.powerW) && noise.powerW>0 && ...
    isfinite(noise.requiredSINR) && noise.requiredSINR>0);
end
