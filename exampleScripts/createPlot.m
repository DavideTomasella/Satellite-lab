addpath("reports\")
load("PALETTE_PLOT1.m")

h = figure(122);
plot(hwp_a(:,1)*pi/180, ...
     hwp_a(:,2),".",MarkerSize=14,Color=MC.blue1)
hold on
plot(0:0.01:pi, ...
    (fitted_data_a(0:0.01:pi)),Color=MC.blue2,LineWidth=1.2)
plot(hwp_b(:,1)*pi/180, ...
     hwp_b(:,2),".",MarkerSize=14,Color=MC.red1)
plot(0:0.01:pi, ...
    (fitted_data_b(0:0.01:pi)),Color=MC.red2,LineWidth=1.2)
hold off
grid on
ylabel("Number of coincidences $N$", ...
       Interpreter="Latex")
xlabel("Waveplate angle $\theta$ [rad]", ...
       Interpreter="Latex")
xlim([0 pi])
xticks([0,pi/6,pi/3,pi/2,2*pi/3,5*pi/6,pi])
xticklabels(["0","\pi/6","\pi/3","\pi/2","2\pi/3","5\pi/6","\pi"])
legend(["","HWP A","","HWP B"], ...
    "Position",[0.65 0.27 0.2 0.1],Interpreter="Latex")
save(h,"reports\tmp")