function [DSq] = DSQUARE(Yin,WINDOW)

% Yin           The data that is to be analysed
% WINDOW        Size of the window for which mu and variance will be
%               calculated before and after the point of interest

DSq = NaN(size(Yin));
for i = (WINDOW+1):(length(Yin)-WINDOW)
    MUb = mean(Yin(i-WINDOW:i-1,1)); %MUb is average from window before point i
    VARb = var(Yin(i-WINDOW:i-1,1));
    MUa = mean(Yin(i+1:i+WINDOW,1)); %MUa is average from window after point i
    VARa = var(Yin(i+1:i+WINDOW,1));
    DSq(i,1) = (MUb-MUa)^2/(VARb+VARa);
end
