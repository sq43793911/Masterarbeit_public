FVec=[0.5,1,2.5,5,7.5,10]; % Kraft in [N]
lVec=[0,16,72,173,274,376]; % Federweg in [mm]
% berechne Ausgleichsgerade
AG = polyfit(lVec,FVec,1);



hold on
% plotte Messpunkte
plot(lVec,FVec,'o')

% plotte Ausgleichsgrade
plot(lVec,AG(2)+AG(1)*lVec,'r-')
box on
grid on
xlabel('l in [mm]')
ylabel('F in [N]')
title('Kraft-Weg-Diagramm Feder')