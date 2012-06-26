% dispNonAff.m 
% ------------
% Makes a graph of the nonaffinity measure.
% This section imports the data and averages.

cd /home/miles/Summer_2012/Summer_Internship/testMiles/
matB = importdata('nonaff_data.txt',',');

matB = sortrows(matB, 3);

p = unique(matB(:, 3));
s = unique(matB(:, 2));
Gamma = zeros(length(p), 1);
j = 1;

for i = 1 : length(p) 
    
    temp = 0.0;
    k = 0;
    
    while matB(j, 3) == p(i) && j < max(size(matB))
        
        temp = temp + matB(j, 4);
        j = j + 1;
        k = k + 1;

    end
    
    Gamma(i) = temp / k;

end

%%

% Plotting the data now.

figure(2)
clf

hold on

ploth = plot(p, Gamma);

pCmarkerx = 0.6484 * ones(100, 1);
pCmarkery = 0:1:99;

plot(pCmarkerx, pCmarkery, '--', 'LineWidth', 3, 'Color', [0.5 0.5 0.5]);

set(gca, ...
    'Box', 'on', ...
    'FontName', 'Nimbus Sans L', ... 
    'FontSize', 16, ...
    'YScale', 'log')

title('Non-affinity measure against bond probability $p$', 'Interpreter', 'latex')

set(ploth, 'Color', [1.0 0.2 0.2], ...
    'LineStyle', 'none', ...
    'LineWidth', 1, ...
    'Marker', 'o', ...
    'MarkerEdgeColor', [1.0, 0.2 0.2], ...h
    'MarkerFaceColor', [0.3 0.3 0.3])

axis([0.45 0.85 0.9 40])

xlabel('Probability of Bond ($p$)', 'Interpreter', 'Latex')
ylabel('Non-affinity Measure', 'Interpreter', 'latex')