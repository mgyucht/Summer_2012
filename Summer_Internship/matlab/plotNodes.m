%% plotNodes.m
% Plot the nodes by their positions

cd /home/miles/Summer_2012/Summer_Internship/integrator/output/
posName = 'position_data';
extension = '.txt';

command = ['ls | grep ', posName, ' | wc -l'];
[dummy, numFiles] = unix(command);
numFiles = str2num(numFiles);

M = moviein(numFiles);

figure(4)

for frame = 1:numFiles

        fileName = strcat(posName, num2str(frame), extension);
    matPos = importdata(fileName,',');

    netSize = matPos.data(1,1);
    strain = matPos.data(1,2);
    shear = strain * netSize * sqrt(3.0)/2.0;

    position = zeros(netSize, netSize, 2);
    springs = zeros(netSize, netSize, 3);

    clf
    hold on
    for i = 0 : netSize - 1
        for j = 0 : netSize - 1
            position(i + 1, j + 1, 1) = matPos.data(2 + (i * netSize + j), 3);
            position(i + 1, j + 1, 2) = matPos.data(2 + (i * netSize + j), 4);
            springs(i+1, j+1, :) = matPos.data(2 + (i * netSize + j), 5:7);
            plot(position(i+1, j+1, 1), position(i+1, j+1, 2), 'ro')
        end
    end

    for i = 1:netSize
        for j = 1:netSize
            if springs(i, j, 1) == 1
                if j == netSize
                    plot([position(i, j, 1), position(i, 1, 1) + netSize], [position(i, j, 2), position(i, 1, 2)], 'b-')
                else
                    plot([position(i, j, 1), position(i, j+1, 1)], [position(i, j, 2), position(i, j+1, 2)], 'b-')
                end
            end
            if springs(i, j, 2) == 1
                if i == netSize
                    plot([position(i, j, 1), position(1, j, 1) + netSize / 2.0 + shear], [position(i, j, 2), position(1, j, 2) + netSize * sqrt(3.0)/2.0], 'b-')
                else
                    plot([position(i, j, 1), position(i+1, j, 1)], [position(i, j, 2), position(i+1, j, 2)], 'b-')
                end
            end
            if springs(i, j, 3) == 1
                if j == 1 && i ~= netSize
                    plot([position(i, j, 1), position(i+1, netSize, 1) - netSize], [position(i, j, 2), position(i+1, netSize, 2)], 'b-')
                elseif i == netSize && j ~= 1
                    plot([position(i, j, 1), position(1, j-1, 1) + netSize / 2.0 + shear], [position(i, j, 2), position(1, j-1, 2) + netSize * sqrt(3.0)/2.0], 'b-')
                elseif i == netSize && j == 1
                    plot([position(i, j, 1), position(1, netSize, 1) - netSize / 2.0 + shear], [position(i, j+1, 2), position(1, netSize-1, 2) + netSize * sqrt(3.0)/2.0], 'b-')
                else
                    plot([position(i, j, 1), position(i+1, j-1, 1)], [position(i, j, 2), position(i+1, j-1, 2)], 'b-')
                end
            end
        end
    end
    axis equal
    title(['Frame ', num2str(frame)])

    M(:, frame) = getframe(gcf);
end

%% Plot one time step

posName = 'position_data';
frame = '';
extension = '.txt';
figure(5)
fileName = strcat(posName, num2str(frame), extension);
matPos = importdata(fileName,',');

netSize = matPos.data(1,1);
strain = matPos.data(1,2);
shear = strain * netSize * sqrt(3.0)/2.0;

position = zeros(netSize, netSize, 2);
springs = zeros(netSize, netSize, 3);

clf
hold on
for i = 0 : netSize - 1
    for j = 0 : netSize - 1
        position(i + 1, j + 1, 1) = matPos.data(2 + (i * netSize + j), 3);
        position(i + 1, j + 1, 2) = matPos.data(2 + (i * netSize + j), 4);
        springs(i+1, j+1, :) = matPos.data(2 + (i * netSize + j), 5:7);
        plot(position(i+1, j+1, 1), position(i+1, j+1, 2), 'ro')
    end
end

for i = 1:netSize
    for j = 1:netSize
        if springs(i, j, 1) == 1
            if j == netSize
                plot([position(i, j, 1), position(i, 1, 1) + netSize], [position(i, j, 2), position(i, 1, 2)], 'b-')
            else
                plot([position(i, j, 1), position(i, j+1, 1)], [position(i, j, 2), position(i, j+1, 2)], 'b-')
            end
        end
        if springs(i, j, 2) == 1
            if i == netSize
                plot([position(i, j, 1), position(1, j, 1) + netSize / 2.0 + shear], [position(i, j, 2), position(1, j, 2) + netSize * sqrt(3.0)/2.0], 'b-')
            else
                plot([position(i, j, 1), position(i+1, j, 1)], [position(i, j, 2), position(i+1, j, 2)], 'b-')
            end
        end
        if springs(i, j, 3) == 1
            if j == 1 && i ~= netSize
                plot([position(i, j, 1), position(i+1, netSize, 1) - netSize], [position(i, j, 2), position(i+1, netSize, 2)], 'b-')
            elseif i == netSize && j ~= 1
                plot([position(i, j, 1), position(1, j-1, 1) + netSize / 2.0 + shear], [position(i, j, 2), position(1, j-1, 2) + netSize * sqrt(3.0)/2.0], 'b-')
            elseif i == netSize && j == 1
                plot([position(i, j, 1), position(1, netSize, 1) - netSize / 2.0 + shear], [position(i, j+1, 2), position(1, netSize-1, 2) + netSize * sqrt(3.0)/2.0], 'b-')
            else
                plot([position(i, j, 1), position(i+1, j-1, 1)], [position(i, j, 2), position(i+1, j-1, 2)], 'b-')
            end
        end
    end
end

title(['Frame ', num2str(frame)])