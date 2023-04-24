function [codeUp, codePhase, dataSize] = upsampleCode(code, f_s, f_c, codePhase)
% Inputs:
%   code    prn
%   f_s     sampling frequency
%   f_c     code frequency
%   samp    number of samples to output
%   cShift  current code chip shift

codePhaseStep = f_c / f_s;

L = length(code);
dataSize = ceil((L - codePhase) / codePhaseStep);

% idx = (1:dataSize).*codePhaseStep + codePhase;
idx = codePhase : codePhaseStep : ((dataSize-1)*codePhaseStep) + codePhase;
idx2 = ceil(idx)+1;

code = [code(end), code, code(1)];

codeUp = code(idx2);
codePhase = idx(dataSize) + codePhaseStep - L;

end