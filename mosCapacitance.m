%%  MATLAB function to calculate capacitances.
%   Date of creation:   09-04-2019
%   Last Modified:      09-04-2019

function [Ciss, Coss, Crss, D_Ciss, D_Coss, D_Crss] = mosCapacitance(mos, vds)

    %%  Calculate Capacitances
    Ciss = mos.ciss(vds)*1e-12;
    Coss = mos.coss(vds)*1e-12;
    Crss = mos.crss(vds)*1e-12;
    
    %%  Calculate Derivatives
    D_Ciss = differentiate(mos.ciss, vds)*1e-12;
    D_Coss = differentiate(mos.coss, vds)*1e-12;
    D_Crss = differentiate(mos.crss, vds)*1e-12;

end
