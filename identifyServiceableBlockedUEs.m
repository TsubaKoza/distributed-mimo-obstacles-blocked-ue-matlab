function classification = identifyServiceableBlockedUEs(result)
%IDENTIFYSERVICEABLEBLOCKEDUES Use actual geometry and retained path lists.
[M,K]=size(result.hasLoS); nlos=false(M,K);
for m=1:M
    for k=1:K
        types=string({result.pathInfo(m,k).paths.type});
        nlos(m,k)=any(types=="reflection" | types=="diffraction");
    end
end
noLoS=all(~result.hasLoS,1);
notOutage=~all(result.outageMask,1);
count=sum(nlos & ~result.outageMask,1);
serviceable=noLoS & count>=2 & notOutage;
classification=struct('allAPBlockedMask',noLoS,'nlosPathMask',nlos, ...
    'nlosAPCount',count,'serviceableMask',serviceable, ...
    'serviceableIndices',find(serviceable),'completeOutageMask',all(result.outageMask,1));
end
