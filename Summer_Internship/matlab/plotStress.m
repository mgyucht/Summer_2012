%% plotStress.m
% Imports the stress vector and plots it against time steps.

cd /home/miles/Summer_2012/Summer_Internship/integrator/output

matStress = importdata('stress_data.txt',',');

figure(3)
clf

hold on

stress_vector = matStress(:, 1);
strain_vector = matStress(:, 2);
time_vector = matStress(:, 3);

plot(time_vector, stress_vector, 'ro')
plot(time_vector, strain_vector, 'b-')

set(gca, ...
    'Box', 'on', ...
    'FontName', 'Nimbus Sans L', ... 
    'FontSize', 16)

title('Shear stress \sigma_{xy} against time steps (\delta t = 0.01)')

xlabel('Time')
ylabel('Stress')
