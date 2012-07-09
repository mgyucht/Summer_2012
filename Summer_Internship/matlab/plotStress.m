%% plotStress.m
% Imports the many stress vectors and calculates shear moduli.

cd /media/sdc1/Summer_2012/Summer_Internship/integrator/output/

clear all;

p_dir = dir;
for i = 3:length(p_dir)
    x = str2num(p_dir(i).name);
    p_val(i-2) = x;
end

for j = 1:length(p_val)
    prob = p_val(j);
    clear w_val;
    clear w_dir;
    cd(p_dir(j + 2).name)
    w_dir = dir;
    for i = 3:length(w_dir)
        w_val(i-2) = str2num(w_dir(i).name);
    end
    
    for i = 1:length(w_val)
        rate = w_val(i);
        cd(w_dir(i + 2).name)
        str_files = dir;
        nFiles = length(str_files) - 2;
        pr_array=zeros(nFiles, 2);
        
        for k = 1:nFiles

            file = strcat('stress_data_', num2str(k), '.txt');
            matStress = importdata(file,',');

        %     figure(3)
        %     clf

        %     hold on

            stress_vector = matStress(:, 1);
            strain_vector = matStress(:, 2);
            time_vector = matStress(:, 3);

        %     plot(time_vector, stress_vector, 'ro')
        %     plot(time_vector, strain_vector, 'b-')
        % 
        %     set(gca, ...
        %         'Box', 'on', ...
        %         'FontName', 'Nimbus Sans L', ... 
        %         'FontSize', 16, ...
        %         'Ylim', [min(strain_vector), max(strain_vector)])
        % 
        %     title('Shear stress \sigma_{xy} against time steps (\delta t = 0.01)')
        % 
        %     xlabel('Time')
        %     ylabel('Stress')

            % Calculate parameters

            a0 = [0.01 rate 0];
            stress_params = lsqcurvefit(@sine_fit, a0, time_vector, stress_vector);
            strain_params = lsqcurvefit(@sine_fit, a0, time_vector, strain_vector);

            stress_amplitude = stress_params(1);
            strain_amplitude = strain_params(1);

            phase_shift = stress_params(3) - strain_params(3);

            G_elastic = stress_amplitude / strain_amplitude * cos(phase_shift);
            G_loss    = stress_amplitude / strain_amplitude * sin(phase_shift);
            pr_array(k, :) = [G_elastic G_loss];
        end
        mean_data = mean(pr_array, 1);
        G_array(i, j, 1:4) = [prob rate mean_data];
        [i j prob rate mean_data]
        cd ..
    end
    cd ..
end

G_array

%% Plot it 3D

figure(6)
clf

surf(G_array1(:, :, 1), G_array1(:, :, 2), G_array1(:, :, 4))
% surf(G_array1(1, :, 1), G_array1(:, 1, 2), G_array1(:, :, 3))

set(gca, 'XScale', 'log', 'YScale', 'log')
axis([min(G_array1(:, 1))*0.9 max(G_array1(:, 1))*1.1 min(G_array1(:, 2:3))*0.9, max(G_array1(:, 2:3))*1.1])

%% Plot it 2D

h9 = figure(9);
clf
hold on
xlabel('Frequency')
ylabel('G''')

for i = 1:8
    x_vec = G_array1(:, i, 2);
    y_vec = G_array1(:, i, 3);
    plot(x_vec, y_vec, 'Color', [(i / 30 + 0.6) 0 0])
    set(gca, 'XScale', 'log', 'Yscale', 'log', 'Ylim', [.001 1])
end

legend('0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8')

saveas(h9, '/home/miles/Summer_2012/Summer_Internship/matlab/G''_vs_f_2d.png')

%%
h8 = figure(8);
clf
hold on
xlabel('Frequency')
ylabel('G''''')

for i = 1:8
    xvec = G_array1(:, i, 2);
    yvec = G_array1(:, i, 4);
    plot(xvec, yvec, 'Color', [(i / 30 + 0.6) 0 0])
    set(gca, 'XScale', 'log', 'Yscale', 'log', 'Ylim', [.01 0.1])
end

legend('0.45', '0.5', '0.55', '0.6', '0.65', '0.7', '0.75', '0.8')
saveas(h8, '/home/miles/Summer_2012/Summer_Internship/matlab/G''''_vs_f_2d.png')
