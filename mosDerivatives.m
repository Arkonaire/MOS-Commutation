%%  MATLAB function to calculate MOSFET current derivatives.
%   Date of creation:   25-03-2019
%   Last Modified:      25-03-2019

function [D_ids, D_id] = mosDerivatives(vgs, vds, D_vgs, D_vds, mos)

    %%  Extract parameters
    Vth = mos(1);
    lambda = mos(2);
    beta = mos(3);
    Is = mos(4);
    Vt = 1.38064852e-23/1.60217662e-19 * (300);
    
    %%  Body diode current derivative calculation
    id = Is*(exp(-vds/Vt) - 1);
    D_id = -((id + Is)/Vt)*D_vds;
    
    %%  MOS current derivative calculation
    if vgs < Vth
        %%  Interdiction region
        D_ids = 0;
    elseif abs(vds) > vgs - Vth
        %%  Saturation region
        D_ids = beta*(vgs - Vth)*(1 + lambda*abs(vds))*sign(vds)*D_vgs + ...
            (beta/2)*((vgs - Vth)^2)*lambda*D_vds;
    else
        %%  Ohmic Region
        D_ids = (beta/2)*(2*(vgs - Vth)*(1 + 2*lambda*abs(vds)))*D_vds - ...
            (beta/2)*abs(vds)*(2 + 3*lambda*abs(vds))*D_vds + ...
            beta*abs(vds)*(1 + lambda*abs(vds))*sign(vds)*D_vgs;
    end

end