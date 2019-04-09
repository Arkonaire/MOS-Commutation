%%  MATLAB script to implement single iteration of forward euler algo to solve stiff ode.
%   Date of creation:   25-03-2019
%   Last Modified:      25-03-2019

function [y, dt] = eulerUpdate(x, mos, abcd, circuit, dt_prev)

    %%  Calculate derivatives
    dx = stateSpaceSRSC(x, mos, abcd, circuit);
    
    %%  Solver parameters
    eIDS = 0.005;
    eId = 0.005;
    eIgs = 0.005;
    eCyx = 0.005;
    alpha = 0.25;
    
    %%  Calculate reference current derivatives
    [D_idsH, D_idH] = mosDerivatives(x(1), x(2), dx(1), dx(2), mos);
    [D_idsL, D_idL] = mosDerivatives(x(3), x(4), dx(3), dx(4), mos);
    D_igsH = -dx(1)/circuit(6);
    D_igsL = -dx(3)/circuit(7);
    
    %%  Calculate reference capacitance and voltage derivatives
    vds = x(2)*[1; 1; 1; 0; 0; 0] + x(4)*[0; 0; 0; 1; 1; 1];
    D_vds = dx(2)*[1 1 1 0 0 0] + dx(4)*[0 0 0 1 1 1];
    D_Cyx = (abcd(:,2).*abcd(:,1).*exp(abcd(:,2).*abs(vds)) + abcd(:,4).*abcd(:,3).*exp(abcd(:,4).*abs(vds)))';
    
    %%  Calculate max variations
    del_ids = 2*circuit(4);
    del_id = circuit(4);
    del_igsH = circuit(1)/circuit(6);
    del_igsL = circuit(2)/circuit(7);
    del_Cyx = (abcd(:,1) + abcd(:,3))';
    
    %%  Calculate reference timestep
    dt_ref = min([eIDS*del_ids/D_idsH, eIDS*del_ids/D_idsL, ...
        eId*del_id/D_idH, eId*del_id/D_idL, ...
        eIgs*del_igsH/D_igsH, eIgs*del_igsL/D_igsL, ...
        eCyx*del_Cyx./(D_vds.*D_Cyx)]);
    
    %%  Calculate max timestep
    CissH = abcd(1,1)*exp(abcd(1,2)*abs(x(2))) + abcd(1,3)*exp(abcd(1,4)*abs(x(2)));
    CissL = abcd(4,1)*exp(abcd(4,2)*abs(x(4))) + abcd(4,3)*exp(abcd(4,4)*abs(x(4)));
    dt_max = 1e-4*min([CissH*circuit(6), CissL*circuit(7)]);
    
    %%  Calculate timestep
    dt = 0.5*(dt_max*(1 + (2/pi)*atan(alpha*(dt_ref - dt_prev)/min(dt_ref, dt_prev))) + dt_prev);
    
    %%  Update data
    dt = 1e-15;
    y = x + dx*dt;

end
