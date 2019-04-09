%%  MATLAB function for modified euler algorithm.
%   Date of creation:   26-03-2019
%   Last Modified:      26-03-2019

function [t, y] = modForwardEuler(mos, ckt, abcd, tspan, ini)

    %%  Initialization
    t0 = tspan(1);
    y0 = ini;
    dt0 = 0;
    t = t0;
    y = y0;
    
    %%  Iterate through tspan
    while(t0 < tspan(2))
        [y0, dt0] = eulerUpdate(y0, mos, abcd, ckt(t0), dt0);
        t0 = t0 + dt0;
        t = [t t0];
        y = [y y0];
    end
    t = t';
    y = y';

end
