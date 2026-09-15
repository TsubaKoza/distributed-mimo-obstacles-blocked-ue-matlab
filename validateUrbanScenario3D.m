function validateUrbanScenario3D(s,cfg)
%VALIDATEURBANSCENARIO3D Independent outdoor/spacing/corridor assertions.
assert(size(s.APpos,1)==cfg.M && size(s.UEpos,1)==cfg.K);
points=[s.APpos;s.UEpos];
assert(all(points(:,1:2)>=cfg.region(1:2,1).','all') && ...
    all(points(:,1:2)<=cfg.region(1:2,2).','all'));
for i=1:numel(s.obstacles)
    o=s.obstacles(i); e=s.buildingEnvelopes(i,:);
    assert(e(1)>=cfg.region(1,1) && e(2)<=cfg.region(1,2) && ...
        e(3)>=cfg.region(2,1) && e(4)<=cfg.region(2,2));
    for q=1:size(points,1)
        local=transformToObstacleLocal(points(q,:),o);
        assert(any(abs(local(1:2))>[o.width,o.depth]/2+cfg.obstacleClearance));
    end
    others=[s.buildingEnvelopes(i+1:end,:);s.corridors];
    for j=1:size(others,1)
        b=others(j,:);
        assert(e(2)<=b(1) || b(2)<=e(1) || e(4)<=b(3) || b(4)<=e(3), ...
            'Building overlaps a building or a building-free corridor.');
    end
end
for m=1:cfg.M
    assert(all(vecnorm(s.APpos(m+1:end,:)-s.APpos(m,:),2,2)>=cfg.minimumAPSeparation));
end
if strcmpi(cfg.reflectionCoefficientMode,"fixed")
    expected=cfg.fixedReflectionMagnitude*exp(1j*deg2rad(cfg.fixedReflectionPhaseDeg));
    assert(all([s.obstacles.reflectionCoefficient]==expected), ...
        'All obstacles must share the configured complex reflection coefficient.');
end
end
