clear
dz = [0.1 0.5 1 2 2.5];
%[0.2 0.4 0.8, 1.6, 3.2];

clf
for i = 1:length(dz)
    [t,z,P] = grid_func(dz(i));
    plot(P(end,:), -z, '--','Linewidth',1.5)
    drawnow
    hold on
end
%%
%legend('0.2','0.4','0.8','1.6','3.2')
legend('0.1','0.5','1','2','2.5')
xlabel('Concentration (X/m^3)')
ylabel('Depth (m)')
title('Grid Size Sensitivity Analysis')

