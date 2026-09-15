function [weights,effectiveAmplitude,effectiveGain] = computeBlockedUEMRT(storedChannel)
%COMPUTEBLOCKEDUEMRT Existing channel g is UE->AP, N_AP x number_of_APs.
% Reciprocal downlink: y = g.'*w*s. hDL=conj(g) gives hDL'*w.
% w=conj(g)/norm(g); effectiveAmplitude=norm(g), real nonnegative.
gain=sum(abs(storedChannel).^2,1);
weights=complex(zeros(size(storedChannel)));
nonzero=gain>0;
weights(:,nonzero)=conj(storedChannel(:,nonzero))./sqrt(gain(nonzero));
effectiveAmplitude=sum(storedChannel.*weights,1); % NON-conjugate transpose
effectiveGain=abs(effectiveAmplitude).^2;
assert(all(abs(effectiveGain(nonzero)-gain(nonzero))<=1e-12*gain(nonzero)));
assert(all(abs(sum(abs(weights(:,nonzero)).^2,1)-1)<1e-12));
end
