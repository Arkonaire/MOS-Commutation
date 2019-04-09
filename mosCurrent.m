%%  MATLAB function to calculate MOSFET currents.
%   Date of creation:   24-03-2019
%   Last Modified:      09-04-2019

function [ids, id] = mosCurrent(mos, vgs, vds)

    %%  Extract parameters [Vth lambda beta Is Rdson]
    Vth = mos.param(1);
    lambda = mos.param(2);
    beta = mos.param(3);
    Is = mos.param(4);
    Vt = 1.38064852e-23/1.60217662e-19 * (300);
    
    %%  Body diode current calculation
    id = Is*(exp(-vds/Vt) - 1);
    
    %%  MOS current calculation
    if vgs < Vth
        %%  Interdiction region
        ids = 0;
    elseif abs(vds) > vgs - Vth
        %%  Saturation region
        ids = (beta/2)*((vgs - Vth)^2)*(1 + lambda*abs(vds))*sign(vds);
    else
        %%  Ohmic Region
        ids = (beta/2)*abs(vds)*(2*(vgs - Vth) - abs(vds))*(1 + lambda*abs(vds))*sign(vds);
    end

end