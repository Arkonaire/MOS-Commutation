%%  MATLAB function to analyse MOSFET commutation.
%   Date of creation:   23-03-2019
%   Last Modified:      07-04-2019

function [tH, yH, tL, yL, loss] = analyzeSRSC(mosH, mosL, circuit)

    %%  MOS Parameters
    RdsonH = mosH.param(5);
    RdsonL = mosL.param(5);
    
    %%  Circuit Parameters [VdrH VdrL Voff Ion Roff RgH RgL Loff LgH LgL freq deadTime simTime]
    VdrH = circuit(1);
    VdrL = circuit(2);
    Voff = circuit(3);
    Ion = circuit(4);
    RgH = circuit(6);
    RgL = circuit(7);
    simTime = circuit(13);
    
    %%  Initial Conditions
    y0_HturnOn = [0; Voff; VdrL; Ion*RdsonL; 0; 0; 0];
    y0_LturnOn = [VdrH; Ion*RdsonH; 0; Voff; Ion; 0; 0];
    
    %%  Solve the odes
    ckt_HturnOn = @(t) [circuit(1)*(t > circuit(12)) 0 circuit(3:10)];
    ckt_LturnOn = @(t) [0 circuit(2)*(t > circuit(12)) circuit(3:10)];
    %[tH, yH] = modForwardEuler(mos, ckt_HturnOn, abcd, [0 simTime], y0_HturnOn);
    %[tL, yL] = modForwardEuler(mos, ckt_LturnOn, abcd, [0 simTime], y0_LturnOn);
    [tH, yH] = ode15s(@(t,y) stateSpaceSRSC(y, mosH, mosL, ckt_HturnOn(t)), [0 simTime], y0_HturnOn, odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    [tL, yL] = ode15s(@(t,y) stateSpaceSRSC(y, mosH, mosL, ckt_LturnOn(t)), [0 simTime], y0_LturnOn, odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    
    %%  Plot VI waveforms
    figure('Name', 'VI Waveforms');
    subplot(3,2,1);
    plot(tH, [yH(:,1) yH(:,3)])
    title('Vgs High side Turn on');
    subplot(3,2,2);
    plot(tL, [yL(:,1) yL(:,3)]);
    title('Vgs Low side Turn on');
    subplot(3,2,3);
    plot(tH, [yH(:,2) yH(:,4)]);
    title('Vds High side Turn on');
    subplot(3,2,4);
    plot(tL, [yL(:,2) yL(:,4)]);
    title('Vds Low side Turn on');
    subplot(3,2,5);
    plot(tH, [yH(:,6) yH(:,7)]);
    title('Ig High side Turn on');
    subplot(3,2,6);
    plot(tL, [yL(:,6) yL(:,7)]);
    title('Ig Low side Turn on');
    
    %%  Format VI plots
    for i = 1:6
        subplot(3,2,i);
        xlabel('Time(s)');
        if i > 4
            ylabel('Current(A)');
        else
            ylabel('Voltage(V)');
        end
        legend('High', 'Low');
        grid on;
    end
    
    %%  Calculate power losses
    P1H = [RgH*yH(:,6).^2 RgL*yH(:,7).^2];
    P1L = [RgH*yL(:,6).^2 RgL*yL(:,7).^2];
    iH = tH;
    iL = tH;
    for i = 1:length(tH)
        [a, b] = mosCurrent(mosH, yH(i,1), yH(i,2));
        iH(i) = a - b;
        [a, b] = mosCurrent(mosL, yH(i,3), yH(i,4));
        iL(i) = a - b;
    end
    P2H = [yH(:,2).*iH yH(:,4).*iL];
    iH = tL;
    iL = tL;
    for i = 1:length(tL)
        [a, b] = mosCurrent(mosH, yL(i,1), yL(i,2));
        iH(i) = a - b;
        [a, b] = mosCurrent(mosL, yL(i,3), yL(i,4));
        iL(i) = a - b;
    end
    P2L = [yL(:,2).*iH yL(:,4).*iL];
    PH = P1H + P2H;
    PL = P1L + P2L;
    
    %%  Plot VI waveforms
    figure('Name', 'Power Waveforms');
    subplot(3,2,1);
    plot(tH, P1H)
    title('Gate Loss High side Turn on');
    subplot(3,2,2);
    plot(tL, P1L);
    title('Gate Loss Low side Turn on');
    subplot(3,2,3);
    plot(tH, P2H);
    title('Drain Loss High side Turn on');
    subplot(3,2,4);
    plot(tL, P2L);
    title('Drain Loss Low side Turn on');
    subplot(3,2,5);
    plot(tH, PH);
    title('Total Loss High side Turn on');
    subplot(3,2,6);
    plot(tL, PL);
    title('Total Loss Low side Turn on');
    
    %%  Format VI plots
    for i = 1:6
        subplot(3,2,i);
        xlabel('Time(s)');
        ylabel('Power(W)');
        legend('High', 'Low');
        grid on;
    end
    
    %%  Evaluate numerical values of losses
    P1H = sum(P1H(2:end, :).*((tH(2:end) - tH(1:end - 1))*[1 1]));
    P1L = sum(P1L(2:end, :).*((tL(2:end) - tL(1:end - 1))*[1 1]));
    P2H = sum(P2H(2:end, :).*((tH(2:end) - tH(1:end - 1))*[1 1]));
    P2L = sum(P2L(2:end, :).*((tL(2:end) - tL(1:end - 1))*[1 1]));
    PH = P1H + P2H;
    PL = P1L + P2L;
    loss = [P1H P1L; P2H P2L; PH PL]*circuit(11);

end
