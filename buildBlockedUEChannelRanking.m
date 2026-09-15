function [ranking,weights,amplitude] = buildBlockedUEChannelRanking(result,k)
%BUILDBLOCKEDUECHANNELRANKING Rank ALL APs; mark candidate eligibility.
g=result.h_true(:,:,k);
[weights,amplitude,effectiveGain]=computeBlockedUEMRT(g);
M=size(g,2); gain=sum(abs(g).^2,1).';
type=repmat("none",M,1); strongest=-Inf(M,1); feasible=false(M,1);
for m=1:M
    paths=result.pathInfo(m,k).paths;
    if isempty(paths), continue; end
    power=abs([paths.complexGain]).^2;
    [value,index]=max(power);
    strongest(m)=10*log10(value); type(m)=paths(index).type;
    types=string({paths.type});
    feasible(m)=~result.hasLoS(m,k) && any(types=="reflection" | types=="diffraction") ...
        && ~result.outageMask(m,k);
end
ranking=table((1:M).',gain,10*log10(gain),effectiveGain.', ...
    10*log10(effectiveGain.'),result.numReflections(:,k), ...
    result.numDiffractions(:,k),type,strongest,feasible,'VariableNames', ...
    {'APIndex','ChannelGainLinear','ChannelGain_dB','EffectiveGainLinear', ...
    'EffectiveGain_dB','NumReflectionPaths','NumDiffractionPaths', ...
    'StrongestPathType','StrongestPathPower_dB','FeasibleNLoSLink'});
ranking=sortrows(ranking,{'EffectiveGainLinear','APIndex'},{'descend','ascend'});
ranking=addvars(ranking,(1:M).','Before',1,'NewVariableNames','Rank');
end
