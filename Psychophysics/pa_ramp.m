function X = pa_ramp(X, N, sampleFreq)
% PA_RAMP(X,N)
%
% Squared-sine onset and squared-cosine offset ramp.
%
% This smooths on- and offset (N samples) of signal X.
%
% PA_RAMP(X, N)
%

% 2013 Marc van Wanrooij
% e-mail: marcvanwanrooij@neural-code.com

if nargin<2
        Fs              = sampleFreq; % Freq (Hz)
        N       = round(5/1000*Fs); % 5 ms
end

X               = X(:);
XLen    = size(X,1);
if (XLen < 2*N)
        error('pandaToolbox:DSP:RampTooLong','Ramp length greater than signal');
else
        RampOn  = sin(0.5*pi*(0:N)/N).^2; % square sine ramp of N samples
        RampOff = fliplr(RampOn);
        head    = 1:(N+1);
        tail    = (XLen-N):XLen;
        for i = 1:size(X,2)
                X(head,i) = RampOn'.*X(head,i);
                tail=round(tail);
                X(tail,i) = RampOff'.*X(tail,i);
        end;
end;