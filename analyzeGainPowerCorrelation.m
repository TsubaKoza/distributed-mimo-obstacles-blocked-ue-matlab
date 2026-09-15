function results = analyzeGainPowerCorrelation(singleTable)
%ANALYZEGAINPOWERCORRELATION No Statistics Toolbox required.
x=singleTable.EffectiveGain_dB; y=singleTable.MinimumPowerFor10Mbps_W;
valid=isfinite(x) & isfinite(y); x=x(valid); y=y(valid); n=numel(x);
pearson=coefficient(x,y); spearman=coefficient(averageRanks(x),averageRanks(y));
results=struct('N',n,'Pearson',pearson,'Spearman',spearman, ...
    'PearsonPValue',NaN,'SpearmanPValue',NaN, ...
    'PValueMethod',"Unavailable: fewer than 3 nonconstant samples", ...
    'Interpretation',"Required power is analytically inverse to effective gain; correlation is not independent empirical evidence. Small N and rejection/target selection limit inference.");
if n>=3 && isfinite(pearson) && isfinite(spearman)
    results.PearsonPValue=pvalue(pearson,n);
    if n<=8
        rx=averageRanks(x); ry=averageRanks(y);
        rx=(rx-mean(rx))/norm(rx-mean(rx));
        ry=(ry-mean(ry))/norm(ry-mean(ry));
        permutations=perms(1:n);
        permuted=reshape(ry(permutations),size(permutations))*rx;
        results.SpearmanPValue=mean(abs(permuted)>=abs(spearman)-1e-12);
        results.PValueMethod="Pearson: two-sided t reference. Spearman: exact two-sided enumeration of all label permutations (N<=8). Both inferential interpretations require exchangeable/independent samples not guaranteed for selected AP links.";
    else
        results.SpearmanPValue=pvalue(spearman,n);
        results.PValueMethod="Pearson: two-sided t reference. Spearman: t approximation (N>8), not exact. Independence and distributional assumptions are not guaranteed for selected AP links.";
    end
end
fprintf('Gain vs required 10 Mbps power: N=%d, Pearson=%g (p=%g), Spearman=%g (p=%g)\n', ...
    n,pearson,results.PearsonPValue,spearman,results.SpearmanPValue);
fprintf('%s\n%s\n',results.PValueMethod,results.Interpretation);
end
function value=coefficient(x,y)
if numel(x)<2, value=NaN; return; end
x=x-mean(x); y=y-mean(y); denominator=norm(x)*norm(y);
if denominator==0, value=NaN; else, value=max(-1,min(1,(x.'*y)/denominator)); end
end
function ranks=averageRanks(x)
[s,order]=sort(x); ranks=zeros(size(x)); i=1;
while i<=numel(x)
    j=i; while j<numel(x) && s(j+1)==s(i), j=j+1; end
    ranks(order(i:j))=(i+j)/2; i=j+1;
end
end
function p=pvalue(r,n)
% Equivalent to a two-sided Student-t tail, using base MATLAB betainc.
p=betainc(max(0,1-r*r),(n-2)/2,0.5);
end
