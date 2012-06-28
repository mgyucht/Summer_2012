% calcShearMod.m
% 0------------0
% 
% Calculates the shear modulus of the spring network given the energies, network
% size, and pBond.

% matA is the original matrix of data.

cd /home/miles/Summer_2012/Summer_Internship/matlab/
matA = importdata('../testMiles/energy_data.txt', ',');
matA = sortrows(matA, 2);

% G holds the calculated shear modulus, and p holds the different probabilities.

s = unique(matA(:, 3));
p = unique(matA(:, 2));
G = zeros(length(p), 1);
j = 1;

for i = 1 : length(p) 
    
    temp = 0.0;
    k = 0;
    [matA(j,2) i]
    
    while matA(j, 2) == p(i) && j < max(size(matA))
        
        temp1 = 2 * matA(j, 1) / (matA(j, 4) * matA(j, 4) * sqrt(3.0) / 2.0 * (matA(j, 3)).^2);
        temp = temp + temp1;
        j = j + 1;
        k = k + 1;

    end
    
    G(i) = temp / k;

end

%% Calculate the complex shear modulus

a0 = [0.005 2 1.571];
[stress_params, resnorm] = lsqcurvefit(@sine_fit, a0, time, stress);
[strain_params, resnorm] = lsqcurvefit(@sine_fit, a0, time, strain);

stress_amplitude = stress_params(1);
strain_amplitude = strain_params(1);

phase_shift = stress_params(3) - strain_params(3);

G_elastic = stress_amplitude / strain_amplitude * cos(phase_shift)
G_loss    = stress_amplitude / strain_amplitude * sin(phase_shift)

%%

% Generate the figure to display the plot.

figh = figure(1);
clf

hold on

ploth = plot(p, G);

pCmarkerx = 0.6484 * ones(1000, 1);
pCmarkery = 0:0.001:0.999;

plot(pCmarkerx, pCmarkery, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

set(gca, ...
    'Box', 'on', ...
    'FontName', 'Nimbus Sans L', ... 
    'FontSize', 16, ...
    'YScale', 'log')

title('Shear Modulus $G$ against bond probability $p$', 'Interpreter', 'latex')

set(ploth, 'Color', [1.0 0.2 0.2], ...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'Marker', 'o', ...
    'MarkerEdgeColor', [1.0, 0.2 0.2], ...h
    'MarkerFaceColor', [0.3 0.3 0.3])

xlabel('Probability of Bond ($p$)', 'Interpreter', 'Latex')
ylabel('Shear Modulus ($G$)', 'Interpreter', 'latex')

