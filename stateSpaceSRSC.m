%%  MATLAB function for SRSC analysis state space equations.
%   Date of Creation:   24-09-2019
%   Last Modified:      09-09-2019

function y = stateSpaceSRSC(x, mosH, mosL, circuit)

    %%  Extract input voltage data
    vgsH = x(1);
    vdsH = x(2);
    vgsL = x(3);
    vdsL = x(4);
    
    %%  Extract external circuit parameters
    vdrH = circuit(1);
    vdrL = circuit(2);
    Voff = circuit(3);
    Ion = circuit(4);
    Roff = circuit(5);
    RgH = circuit(6);
    RgL = circuit(7);
    
    %%  Extract inductances and currents if specified
    if length(x) > 4
        idsH = x(5);
        igH = x(6);
        igL = x(7);
        Loff = circuit(8);
        LgL = circuit(9);
        LgH = circuit(10);
    else
        idsH = (Voff - vdsH - vdsL)/Roff;
        igH = (vdrH - vgsH)/RgH;
        igL = (vdrL - vgsL)/RgL;
    end
    
    %%  Calculate capacitances and their derivatives
    [CissH, CossH, CrssH, D_CissH, ~, D_CrssH] = mosCapacitance(mosH, vdsH);
    [CissL, CossL, CrssL, D_CissL, ~, D_CrssL] = mosCapacitance(mosL, vdsL);
    
    %%  Generate gamma coefficients
    gamma_gsH = CissH - CrssH;
    gamma_gdH = CrssH;
    gamma_dsH = 0;
    gamma_gsL = CissL - CrssL;
    gamma_gdL = CrssL;
    gamma_dsL = 0;
    
    %%  Generate delta coefficients
    delta_gsH = (D_CissH - D_CrssH)*sign(vdsH)*vgsH;
    delta_gdH = sign(vdsH)*(vgsH - vdsH)*D_CrssH - CrssH;
    delta_dsH = CossH - CrssH;
    delta_gsL = (D_CissL - D_CrssL)*sign(vdsL)*vgsL;
    delta_gdL = sign(vdsL)*(vgsL - vdsL)*D_CrssL - CrssL;
    delta_dsL = CossL - CrssL;
    
    %%  Generate final set of coefficients
    delta_ggH = delta_gsH + delta_gdH;
    delta_dgH = delta_dsH - delta_gdH;
    gamma_ggH = gamma_gsH + gamma_gdH;
    gamma_dgH = gamma_dsH - gamma_gdH;
    delta_ggL = delta_gsL + delta_gdL;
    delta_dgL = delta_dsL - delta_gdL;
    gamma_ggL = gamma_gsL + gamma_gdL;
    gamma_dgL = gamma_dsL - gamma_gdL;
    lambda_H = 1/(delta_dgH*gamma_ggH - delta_ggH*gamma_dgH);
    lambda_L = 1/(delta_dgL*gamma_ggL - delta_ggL*gamma_dgL);
    
    %%  Calculate currents
    [iDSH, idH] = mosCurrent(mosH, vgsH, vdsH);
    [iDSL, idL] = mosCurrent(mosL, vgsL, vdsL);
    
    %%  Initialize output
    y = zeros(length(x),1);
    
    %%  Calculate voltage derivatives
    y(1) = +lambda_H*(delta_dgH*igH + delta_ggH*(iDSH - idH - idsH));
    y(2) = -lambda_H*(gamma_dgH*igH + gamma_ggH*(iDSH - idH - idsH));
    y(3) = +lambda_L*(delta_dgL*igL + delta_ggL*(iDSL - idL - idsH + Ion));
    y(4) = -lambda_L*(gamma_dgL*igL + gamma_ggL*(iDSL - idL - idsH + Ion));
    
    %%  Calculate current derivatives if necessary
    if length(y) > 4
        y(5) = (Voff - vdsH - vdsL - idsH*Roff)/Loff;
        y(6) = (vdrH - vgsH - igH*RgH)/LgH;
        y(7) = (vdrL - vgsL - igL*RgL)/LgL;
    end

end
