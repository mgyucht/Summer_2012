%% plotStress.m
% Imports the stress vector and plots it against time steps.

cd /home/miles/Summer_2012/Summer_Internship/integrator/
matStress = importdata('stress_data.txt',',');

figure(3)
clf

hold on

x = 0:999;
y = 0.1 * sin(10 * x * 0.01);

plot(matStress(:, 2), matStress(:, 1), 'ro')
plot(x, y)

set(gca, ...
    'Box', 'on', ...
    'FontName', 'Nimbus Sans L', ... 
    'FontSize', 16)

title('Shear stress \sigma_{xy} against time steps (\delta t = 0.01)')

xlabel('Time Steps')
ylabel('Stress')