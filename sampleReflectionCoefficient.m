function coefficient = sampleReflectionCoefficient(cfg)
%SAMPLEREFLECTIONCOEFFICIENT Fixed material baseline or legacy random draws.
if strcmpi(cfg.reflectionCoefficientMode,"fixed")
    coefficient = cfg.fixedReflectionMagnitude*exp(1j*deg2rad(cfg.fixedReflectionPhaseDeg));
else
    bounds = cfg.obstacleReflectionMagnitudeRange;
    magnitude = bounds(1)+diff(bounds)*rand;
    coefficient = magnitude*exp(1j*2*pi*rand);
end
end
