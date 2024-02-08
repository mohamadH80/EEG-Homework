%% Init
clc;clear;close all;

%% A(الف)
% Define parameters
ModelParams.R = [8 8.5 9.2] ;
ModelParams.Sigma = [3.3e-3 8.25e-5 3.3e-3];
ModelParams.Lambda = [.5979 .2037 .0237];
ModelParams.Mu = [.6342 .9364 1.0362];

Resolution = 1;
[LocMat,GainMat] = ForwardModel_3shell(Resolution, ModelParams);

% Plotting the dipoles in 3D space
figure;
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'filled');
title('Dipoles in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');

% Save the Gain Matrix (Lead-Field Matrix)
save('GainMatrix.mat', 'GainMat');


%% B(ب)
close;
load('ElecPosXYZ.mat');

% Extract electrode positions and names
ElectrodePos = [];
ElectrodeNames = {};
for i = 1:length(ElecPos)
    ElectrodePos(i, :) = ModelParams.R(3) * ElecPos{i}.XYZ; 
    ElectrodeNames{i} = ElecPos{i}.Name;
end

figure;
% Plot dipoles
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'filled');
title('Dipoles and Electrodes in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
hold on;
% Plot electrodes
scatter3(ElectrodePos(:, 1), ElectrodePos(:, 2), ElectrodePos(:, 3), 'filled', 'r');

% Annotate electrodes
for i = 1:length(ElectrodeNames)
    text(ElectrodePos(i, 1), ElectrodePos(i, 2), ElectrodePos(i, 3), ElectrodeNames{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold off;
%% C(پ)
% Identify surface dipoles
surfaceDipolesIndices = find(vecnorm(LocMat) >= ModelParams.R(3)*0.85);
numSurfaceDipoles = length(surfaceDipolesIndices);

% Select a random surface dipole
randomSurfaceIndex = randi(numSurfaceDipoles); 
selectedDipoleIndex = surfaceDipolesIndices(randomSurfaceIndex);
selectedDipole = LocMat(:, selectedDipoleIndex);

% Compute radial direction (normalize vector)
radialDirection = selectedDipole / norm(selectedDipole);

% Plot the selected dipole distinguishably
hold on;
scatter3(selectedDipole(1), selectedDipole(2), selectedDipole(3), 100, 'm', 'filled');
quiver3(0, 0, 0, radialDirection(1), radialDirection(2), radialDirection(3), norm(selectedDipole), 'm', 'LineWidth', 2);
hold off;


%% D(ت)
clc;close;
load('Interictal.mat'); 
selectedRow = Interictal(1, :);
numDipoles = length(LocMat);
numElctrodes = size(GainMat,1);
Q = zeros(3*numDipoles, length(Interictal));
Q(3*selectedDipoleIndex-2,:) = selectedRow*radialDirection(1);
Q(3*selectedDipoleIndex-1,:) = selectedRow*radialDirection(2);
Q(3*selectedDipoleIndex,  :) = selectedRow*radialDirection(3);
electrodePotentials = GainMat * Q;

figure;
timePoints = 1:size(electrodePotentials, 2);

% Iterate over each electrode
for i = 1:numElctrodes
    subplot(7, 3, i);
    plot(electrodePotentials(i, :));
    title(ElectrodeNames{i});    
    xlim([timePoints(1), timePoints(end)]);
end

%% E(ث)
close; clc;
% Find peaks in the selectedRow
[peaks, locations] = findpeaks(selectedRow,"MinPeakHeight",5);

% Define the window size
windowSize = 7; % Total points in the window
halfWindow = floor(windowSize / 2);

% Initialize the vector for average potentials
averagePotentials = zeros(1, 21);

% Calculate the average potential for each electrode in the windows
for i = 1:numElctrodes 
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    % Calculate the average potential
    averagePotentials(i) = electrode_peak_potential/(length(locations)*windowSize);
end

% Normalize the potentials for color mapping
normalizedPotentials = (averagePotentials - min(averagePotentials)) / (max(averagePotentials) - min(averagePotentials));

% Generate colors based on normalized potentials
colormap jet;  % You can also try 'hot' or other colormaps
colors = colormap(jet(length(normalizedPotentials)));

% Plot dipoles
scatter3(ElectrodePos(:, 1), ElectrodePos(:, 2), ElectrodePos(:, 3), 100, colors, 'filled');
title('Electrode Potentials in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
hold on;

% Annotate electrodes
for i = 1:length(ElectrodeNames)
    text(ElectrodePos(i, 1), ElectrodePos(i, 2), ElectrodePos(i, 3), ElectrodeNames{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
colorbar;
hold off;

%% F(ج)
Display_Potential_3D(ModelParams.R(3), normalizedPotentials)

%% G(چ)
clc;
alpha = 1;
Q_MNE = GainMat'*inv(GainMat*GainMat'+alpha*eye(numElctrodes))*averagePotentials';

Omega = zeros(numDipoles, numDipoles);
for i = 1:numDipoles
    temp = 0;
    for j = 1:numElctrodes
        temp = temp + GainMat(j, 3*i-2:3*i)*GainMat(j, 3*i-2:3*i)';
    end
    Omega(i, i) = sqrt(temp);
end
W = kron(Omega, eye(3));
Wm = inv(W'*W);
alpha_w = 0.01;
Q_WMNE = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_w*eye(numElctrodes))*averagePotentials';

%% H(ح)
clc;
Q_MNE_norms = zeros(numDipoles,1);
Q_WMNE_norms = zeros(numDipoles,1);

for i=1:numDipoles
    Q_MNE_norms(i) = norm(Q_MNE(3*i-2:3*i));
    Q_WMNE_norms(i) = norm(Q_WMNE(3*i-2:3*i));
end
[mne_max,mne_max_index] = max(Q_MNE_norms);
[wmne_max,wmne_max_index] = max(Q_WMNE_norms);
% MNE method maxDipole specifications
mne_MS = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
mne_MS(1) = mne_max_index;
mne_MS(2) = mne_max;
mne_MS(3:5) = LocMat(:,mne_max_index);
mne_MS(6:8) = Q_MNE(3*mne_max_index-2:3*mne_max_index);

wmne_MS = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
wmne_MS(1) = wmne_max_index;
wmne_MS(2) = wmne_max;
wmne_MS(3:5) = LocMat(:,wmne_max_index);
wmne_MS(6:8) = Q_MNE(3*wmne_max_index-2:3*wmne_max_index);

% Printing specifications
disp('MNE Maximum Dipole Specifications:');
disp('index:')
disp(mne_MS(1));
disp('norm:')
disp(mne_MS(2));
disp('location:')
disp(mne_MS(3:5));
disp('direction:')
disp(mne_MS(6:8));
disp('wMNE Maximum Dipole Specifications:');
disp('index:')
disp(wmne_MS(1));
disp('norm:')
disp(wmne_MS(2));
disp('location:')
disp(wmne_MS(3:5));
disp('direction:')
disp(wmne_MS(6:8));

%% I(خ)
clc;
real_loc = selectedDipole;
mne_location_error = norm(real_loc-mne_MS(3:5));
mne_direction_error = acosd(mne_MS(6:8)'*real_loc/norm(real_loc));
wmne_location_error = norm(real_loc-wmne_MS(3:5));
wmne_direction_error = acosd(wmne_MS(6:8)'*real_loc/norm(real_loc));
% Print
disp('MNE Maximum Dipole Error:');
disp('distance error:')
disp(mne_location_error)
disp('direction error:')
disp(mne_direction_error)
disp('wMNE Maximum Dipole Error:');
disp('distance error:')
disp(wmne_location_error)
disp('direction error:')
disp(wmne_direction_error)

%% J(د)--B
clc;close all;
figure;
% Plot dipoles
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'filled');
title('Dipoles and Electrodes in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
hold on;
% Plot electrodes
scatter3(ElectrodePos(:, 1), ElectrodePos(:, 2), ElectrodePos(:, 3), 'filled', 'r');

% Annotate electrodes
for i = 1:length(ElectrodeNames)
    text(ElectrodePos(i, 1), ElectrodePos(i, 2), ElectrodePos(i, 3), ElectrodeNames{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end

hold off;

%% J(د)--C
% Identify surface dipoles
surfaceDipolesIndices = find(vecnorm(LocMat) < ModelParams.R(3)*0.5);
numSurfaceDipoles = length(surfaceDipolesIndices);

% Select a random surface dipole
randomSurfaceIndex = randi(numSurfaceDipoles); 
selectedDipoleIndex = surfaceDipolesIndices(randomSurfaceIndex);
selectedDipole = LocMat(:, selectedDipoleIndex);

% Compute radial direction (normalize vector)
radialDirection = selectedDipole / norm(selectedDipole);

% Plot the selected dipole distinguishably
hold on;
scatter3(selectedDipole(1), selectedDipole(2), selectedDipole(3), 100, 'm', 'filled');
quiver3(0, 0, 0, radialDirection(1), radialDirection(2), radialDirection(3), norm(selectedDipole), 'm', 'LineWidth', 2);
hold off;

%% J(د)--D
clc;close;
load('Interictal.mat'); 
selectedRow = Interictal(1, :);
Q = zeros(3*numDipoles, length(Interictal));
randomIndex = selectedDipoleIndex;
Q(3*randomIndex-2,:) = selectedRow*radialDirection(1);
Q(3*randomIndex-1,:) = selectedRow*radialDirection(2);
Q(3*randomIndex,  :) = selectedRow*radialDirection(3);
electrodePotentials = GainMat * Q;

figure;
timePoints = 1:size(electrodePotentials, 2); % or use actual time values

% Iterate over each electrode
for i = 1:numElctrodes
    % Create a subplot for each electrode
    subplot(7, 3, i);
    plot(electrodePotentials(i, :));
    title(ElectrodeNames{i});    
    xlim([timePoints(1), timePoints(end)]);
end
%% J(د)--E
close; clc;
% Find peaks in the selectedRow
[peaks, locations] = findpeaks(selectedRow,"MinPeakHeight",5);

% Define the window size
windowSize = 7; % Total points in the window
halfWindow = floor(windowSize / 2);

% Initialize the vector for average potentials
averagePotentials = zeros(1, numElctrodes);

% Calculate the average potential for each electrode in the windows
for i = 1:21  % For each electrode
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    % Calculate the average potential
    averagePotentials(i) = electrode_peak_potential/(length(locations)*windowSize);
end

% Normalize the potentials for color mapping
normalizedPotentials = (averagePotentials - min(averagePotentials)) / (max(averagePotentials) - min(averagePotentials));

% Generate colors based on normalized potentials
colormap jet;  % You can also try 'hot' or other colormaps
colors = colormap(jet(length(normalizedPotentials)));

% Plot dipoles
scatter3(ElectrodePos(:, 1), ElectrodePos(:, 2), ElectrodePos(:, 3), 100, colors, 'filled');
title('Electrode Potentials in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
hold on;

% Annotate electrodes
for i = 1:length(ElectrodeNames)
    text(ElectrodePos(i, 1), ElectrodePos(i, 2), ElectrodePos(i, 3), ElectrodeNames{i}, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
end
colorbar;
hold off;

%% J(د)--F
Display_Potential_3D(ModelParams.R(3), normalizedPotentials)

%% J(د)--G
clc;
alpha = 1;
Q_MNE = GainMat'*inv(GainMat*GainMat'+eye(21))*averagePotentials';

Omega = zeros(numDipoles, numDipoles);
for i = 1:numDipoles
    temp = 0;
    for j = 1:numElctrodes
        temp = temp + GainMat(j, 3*i-2:3*i)*GainMat(j, 3*i-2:3*i)';
    end
    Omega(i, i) = sqrt(temp);
end
W = kron(Omega, eye(3));
Wm = inv(W'*W);
alpha_w = 0.01;
Q_WMNE = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_w*eye(21))*averagePotentials';

%% J(د)--H
clc;
Q_MNE_norms = zeros(numDipoles,1);
Q_WMNE_norms = zeros(numDipoles,1);

for i=1:numDipoles
    Q_MNE_norms(i) = norm(Q_MNE(3*i-2:3*i));
    Q_WMNE_norms(i) = norm(Q_WMNE(3*i-2:3*i));
end
[mne_max,mne_max_index] = max(Q_MNE_norms);
[wmne_max,wmne_max_index] = max(Q_WMNE_norms);
% MNE method maxDipole specifications
mne_MS = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
mne_MS(1) = mne_max_index;
mne_MS(2) = mne_max;
mne_MS(3:5) = LocMat(:,mne_max_index);
mne_MS(6:8) = Q_MNE(3*mne_max_index-2:3*mne_max_index);

wmne_MS = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
wmne_MS(1) = wmne_max_index;
wmne_MS(2) = wmne_max;
wmne_MS(3:5) = LocMat(:,wmne_max_index);
wmne_MS(6:8) = Q_MNE(3*wmne_max_index-2:3*wmne_max_index);

% Printing specifications
disp('MNE Maximum Dipole Specifications:');
disp('index:')
disp(mne_MS(1));
disp('norm:')
disp(mne_MS(2));
disp('location:')
disp(mne_MS(3:5));
disp('direction:')
disp(mne_MS(6:8));
disp('wMNE Maximum Dipole Specifications:');
disp('index:')
disp(wmne_MS(1));
disp('norm:')
disp(wmne_MS(2));
disp('location:')
disp(wmne_MS(3:5));
disp('direction:')
disp(wmne_MS(6:8));

%% J(د)--I
clc;
real_loc = LocMat(:,randomIndex);
mne_location_error = norm(LocMat(:,randomIndex)-mne_MS(3:5));
mne_direction_error = acosd(-mne_MS(6:8)'*real_loc/norm(real_loc));
wmne_location_error = norm(LocMat(:,randomIndex)-wmne_MS(3:5));
wmne_direction_error = acosd(-wmne_MS(6:8)'*real_loc/norm(real_loc));
% Print
disp('MNE Maximum Dipole Error:');
disp('distance error:')
disp(mne_location_error)
disp('direction error:')
disp(mne_direction_error)
disp('wMNE Maximum Dipole Error:');
disp('distance error:')
disp(wmne_location_error)
disp('direction error:')
disp(wmne_direction_error)

%% K(ذ)--G
clc;
surfaceDipolesIndices = find(vecnorm(LocMat) < ModelParams.R(3)*0.5);
deepDipolesIndices = find(vecnorm(LocMat) > ModelParams.R(3)*0.85);
numSurfaceDipoles = length(surfaceDipolesIndices);
numDeepDipoles = length(deepDipolesIndices);

randomSurfaceIndex_loreta = randi(numSurfaceDipoles); 
randomDeepIndex_loreta = randi(numDeepDipoles);
selectedDipoleIndex_loreta_surface = surfaceDipolesIndices(randomSurfaceIndex_loreta);
selectedDipoleIndex_loreta_deep = deepDipolesIndices(randomDeepIndex_loreta);

selectedDipole_loreta_surface = LocMat(:, selectedDipoleIndex_loreta_surface);
selectedDipole_loreta_deep = LocMat(:, selectedDipoleIndex_loreta_deep);

radialDirection_surface = selectedDipole_loreta_surface / norm(selectedDipole_loreta_surface);
radialDirection_deep = selectedDipole_loreta_deep / norm(selectedDipole_loreta_deep);

Q_surface = zeros(3*numDipoles, length(Interictal));
Q_surface(3*selectedDipoleIndex_loreta_surface-2,:) = selectedRow*radialDirection_surface(1);
Q_surface(3*selectedDipoleIndex_loreta_surface-1,:) = selectedRow*radialDirection_surface(2);
Q_surface(3*selectedDipoleIndex_loreta_surface,  :) = selectedRow*radialDirection_surface(3);
electrodePotentials_surface = GainMat * Q_surface;
Q_deep = zeros(3*numDipoles, length(Interictal));
Q_deep(3*selectedDipoleIndex_loreta_deep-2,:) = selectedRow*radialDirection_deep(1);
Q_deep(3*selectedDipoleIndex_loreta_deep-1,:) = selectedRow*radialDirection_deep(2);
Q_deep(3*selectedDipoleIndex_loreta_deep,  :) = selectedRow*radialDirection_deep(3);
electrodePotentials_deep = GainMat * Q_deep;

windowSize = 7;
halfWindow = floor(windowSize / 2);
averagePotentials_surface = zeros(1, 21);
for i = 1:numElctrodes 
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials_surface(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    averagePotentials_surface(i) = electrode_peak_potential/(length(locations)*windowSize);
end
averagePotentials_deep = zeros(1, numElctrodes);
for i = 1:numElctrodes 
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials_surface(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    averagePotentials_deep(i) = electrode_peak_potential/(length(locations)*windowSize);
end

%--------------------------------------
A1 = zeros(numDipoles, numDipoles);
minnorm = 1000;
for i=1:numDipoles
    for j=1:numDipoles
        if norm(LocMat(:,i)-LocMat(:,j)) == 1
            A1(i,j) = 1/6;
        end
    end
end
onep = ones(numDipoles,1);
A0 = inv(diag(A1*onep)+eye(numDipoles))*A1;
A = kron(A0, eye(3));
B = 6*(A-eye(3*numDipoles));
Omega = zeros(numDipoles, numDipoles);
for i = 1:numDipoles
    temp = 0;
    for j = 1:numElctrodes
        temp = temp + GainMat(j, 3*i-2:3*i)*GainMat(j, 3*i-2:3*i)';
    end
    Omega(i, i) = sqrt(temp);
end
Ww = kron(Omega, eye(3));
W = Ww * B' * B * Ww;
Wm = inv(W'*W);
alpha_loreta = 0.01;
Q_loreta_surface = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_loreta*eye(21))*averagePotentials_surface';
Q_loreta_deep = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_loreta*eye(21))*averagePotentials_deep';

%--------------------------------------------------------------------------
CovarianceMatrix = electrodePotentials_surface*electrodePotentials_surface'/length(electrodePotentials_surface);
Cin = inv(CovarianceMatrix);
alpha_sloreta = 0.01;
Q_sloreta_surface = inv(GainMat'*Cin*GainMat+alpha_sloreta*eye(3*numDipoles))*GainMat'*Cin*averagePotentials_surface';
Q_sloreta_deep = inv(GainMat'*Cin*GainMat+alpha_sloreta*eye(3*numDipoles))*GainMat'*Cin*averagePotentials_deep';


%% K(ذ)--H
clc;
Q_loreta_norms_surface = zeros(numDipoles,1);
Q_loreta_norms_deep = zeros(numDipoles,1);
Q_sloreta_norms_surface = zeros(numDipoles,1);
Q_sloreta_norms_deep = zeros(numDipoles,1);

for i=1:numDipoles
    Q_loreta_norms_surface(i) = norm(Q_loreta_surface(3*i-2:3*i));
    Q_loreta_norms_deep(i) = norm(Q_loreta_deep(3*i-2:3*i));
    Q_sloreta_norms_surface(i) = norm(Q_sloreta_surface(3*i-2:3*i));
    Q_sloreta_norms_surface(i) = norm(Q_sloreta_surface(3*i-2:3*i));
end
[loreta_max_surface,loreta_max_index_surface] = max(Q_loreta_norms_surface);
[loreta_max_deep,loreta_max_index_deep] = max(Q_loreta_norms_deep);
[sloreta_max_surface,sloreta_max_index_surface] = max(Q_sloreta_norms_surface);
[sloreta_max_deep,sloreta_max_index_deep] = max(Q_sloreta_norms_deep);


loreta_MS_surface = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
loreta_MS_surface(1) = loreta_max_index_surface;
loreta_MS_surface(2) = loreta_max_surface;
loreta_MS_surface(3:5) = LocMat(:,loreta_max_index_surface);
loreta_MS_surface(6:8) = Q_loreta_surface(3*loreta_max_index_surface-2:3*loreta_max_index_surface);

loreta_MS_deep = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
loreta_MS_deep(1) = loreta_max_index_deep;
loreta_MS_deep(2) = loreta_max_deep;
loreta_MS_deep(3:5) = LocMat(:,loreta_max_index_deep);
loreta_MS_deep(6:8) = Q_loreta_deep(3*loreta_max_index_deep-2:3*loreta_max_index_deep);


sloreta_MS_surface = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
sloreta_MS_surface(1) = sloreta_max_index_surface;
sloreta_MS_surface(2) = sloreta_max_surface;
sloreta_MS_surface(3:5) = LocMat(:,sloreta_max_index_surface);
sloreta_MS_surface(6:8) = Q_sloreta_surface(3*sloreta_max_index_surface-2:3*sloreta_max_index_surface);

sloreta_MS_deep = zeros(8,1); % index, norm, location(x,y,z), direction(dx,dy,dz)
sloreta_MS_deep(1) = sloreta_max_index_deep;
sloreta_MS_deep(2) = sloreta_max_deep;
sloreta_MS_deep(3:5) = LocMat(:,sloreta_max_index_deep);
sloreta_MS_deep(6:8) = Q_sloreta_deep(3*sloreta_max_index_deep-2:3*sloreta_max_index_deep);


% Printing specifications
disp('Loreta Maximum Dipole Specifications(Surface Dipole):');
disp('index:')
disp(loreta_MS_surface(1));
disp('norm:')
disp(loreta_MS_surface(2));
disp('location:')
disp(loreta_MS_surface(3:5));
disp('direction:')
disp(loreta_MS_surface(6:8));
% Printing specifications
disp('Loreta Maximum Dipole Specifications(Deep Dipole):');
disp('index:')
disp(loreta_MS_deep(1));
disp('norm:')
disp(loreta_MS_deep(2));
disp('location:')
disp(loreta_MS_deep(3:5));
disp('direction:')
disp(loreta_MS_deep(6:8));


% Printing specifications
disp('sLoreta Maximum Dipole Specifications(Surface Dipole):');
disp('index:')
disp(sloreta_MS_surface(1));
disp('norm:')
disp(sloreta_MS_surface(2));
disp('location:')
disp(sloreta_MS_surface(3:5));
disp('direction:')
disp(sloreta_MS_surface(6:8));
% Printing specifications
disp('sLoreta Maximum Dipole Specifications(Deep Dipole):');
disp('index:')
disp(sloreta_MS_deep(1));
disp('norm:')
disp(sloreta_MS_deep(2));
disp('location:')
disp(sloreta_MS_deep(3:5));
disp('direction:')
disp(sloreta_MS_deep(6:8));

%% K(ذ)--I
clc;
real_loc = LocMat(:,selectedDipoleIndex_loreta_surface);
loreta_location_error = norm(real_loc-loreta_MS_surface(3:5));
loreta_direction_error = acosd(loreta_MS_surface(6:8)'*real_loc/norm(real_loc));
% Print
disp('Loreta Maximum Dipole Error(Surface):');
disp('distance error:')
disp(loreta_location_error)
disp('direction error:')
disp(loreta_direction_error)


real_loc = LocMat(:,selectedDipoleIndex_loreta_deep);
loreta_location_error = norm(real_loc-loreta_MS_deep(3:5));
loreta_direction_error = acosd(loreta_MS_deep(6:8)'*real_loc/norm(real_loc));
% Print
disp('Loreta Maximum Dipole Error(Deep):');
disp('distance error:')
disp(loreta_location_error)
disp('direction error:')
disp(loreta_direction_error)
%-----------------------------------
real_loc = LocMat(:,selectedDipoleIndex_loreta_surface);
sloreta_location_error = norm(real_loc-sloreta_MS_surface(3:5));
sloreta_direction_error = acosd(sloreta_MS_surface(6:8)'*real_loc/norm(real_loc));
% Print
disp('sLoreta Maximum Dipole Error(Surface):');
disp('distance error:')
disp(sloreta_location_error)
disp('direction error:')
disp(sloreta_direction_error)


real_loc = LocMat(:,selectedDipoleIndex_loreta_deep);
sloreta_location_error = norm(real_loc-sloreta_MS_deep(3:5));
sloreta_direction_error = acosd(sloreta_MS_deep(6:8)'*real_loc/norm(real_loc));
% Print
disp('sLoreta Maximum Dipole Error(Deep):');
disp('distance error:')
disp(sloreta_location_error)
disp('direction error:')
disp(sloreta_direction_error)



%% L(ر)
clc;close all;
% Identify surface dipoles
surfaceDipolesIndices = find(vecnorm(LocMat) < ModelParams.R(3)*0.6);
numSurfaceDipoles = length(surfaceDipolesIndices);

% Select a random surface dipole
randomSurfaceIndex = randi(numSurfaceDipoles); 
selectedDipoleIndex = surfaceDipolesIndices(randomSurfaceIndex);
selectedDipole = LocMat(:, selectedDipoleIndex);

% Compute radial direction (normalize vector)
radialDirection = selectedDipole / norm(selectedDipole);


selectedRow = Interictal(1, :);
Q = zeros(3*numDipoles, length(Interictal));
Q(3*selectedDipoleIndex-2,:) = selectedRow*radialDirection(1);
Q(3*selectedDipoleIndex-1,:) = selectedRow*radialDirection(2);
Q(3*selectedDipoleIndex,  :) = selectedRow*radialDirection(3);
electrodePotentials = GainMat * Q;

% MUSIC Algorithm for EEG/MEG Source Localization

% GainMat should be of size [Number_of_Sensors, Number_of_Sources]
% Data should be of size [Number_of_Sensors, Number_of_Time_Points]

% Step 1: Compute the data covariance matrix
R = cov(electrodePotentials');

% Step 2: Perform Singular Value Decomposition (SVD)
[U, S, ~] = svd(R);

% Step 3: Determine the number of sources (p)
% This can be a predefined number or determined from the data
p = 1; % Example: assuming 2 sources

% Step 4: Form the signal and noise subspace
Us = U(:, 1:p); 
Un = U(:, p+1:end);

% Step 5: MUSIC Scan
numDipoles = size(GainMat, 2) / 3; % Assuming 3 components per dipole
music_spectrum = zeros(numDipoles, 1);

for i = 1:numDipoles
    lead_field = GainMat(:, 3*i-2:3*i);
    music_spectrum(i) = 1 / norm(lead_field' * Un, 'fro')^2;
end

% Step 6: Identify Peaks in the MUSIC Spectrum
% This identifies the most likely source locations
[peaks, locs] = findpeaks(music_spectrum,'MinPeakDistance',100);
mean_peaks = mean(peaks);
std_peaks = std(peaks);
threshold = mean_peaks + std_peaks; % One standard deviation above the mean
source_locs = locs(peaks > threshold); % Define a suitable threshold

% Step 7: Plot Results
plot(music_spectrum);
hold on;
plot(locs, peaks, 'r*'); % Marking the peaks
title('MUSIC Spectrum');
xlabel('Dipole Index');
ylabel('Spectral Power');
hold off;

% source_locs contains the indices of the estimated source locations
estimated_locations = LocMat(:, locs);  % Retrieves the 3D coordinates



% Assuming true_locs and estimated_locs are Nx3 matrices (N sources, x/y/z coordinates)
% Calculate Euclidean distance for each source
location_errors = sqrt(sum((selectedDipole - estimated_locations).^2, 2))


%% M(ز)
selected_dipoles_indices = [714, 715, 728, 729, 880, 881, 893, 894, 895, 896, 906, 907, 908, 909, 919, 920, 1050, 1051, 1062, 1063];

figure;
% Plot dipoles
scatter3(LocMat(1, :), LocMat(2, :), LocMat(3, :), 'filled');
title('Dipoles and Electrodes in 3D Space');
xlabel('X (cm)');
ylabel('Y (cm)');
zlabel('Z (cm)');
hold on;

% Plot the selected dipole distinguishably
for i=1:length(selected_dipoles_indices)
    selectedDipole = LocMat(:,selected_dipoles_indices(i));
    radialDirection = selectedDipole / norm(selectedDipole);
    scatter3(selectedDipole(1), selectedDipole(2), selectedDipole(3), 100, 'm', 'filled');
    quiver3(0, 0, 0, radialDirection(1), radialDirection(2), radialDirection(3), norm(selectedDipole), 'm', 'LineWidth', 2);
    hold on;
end

%% N(ژ)
for i=1:length(selected_dipoles_indices)
    selectedRow = Interictal(i, :);
    selectedDipoleIndex = selected_dipoles_indices(i);
    radialDirection = LocMat(:, selectedDipoleIndex)/norm(LocMat(:, selectedDipoleIndex));
    Q = zeros(3*numDipoles, length(Interictal));
    Q(3*selectedDipoleIndex-2,:) = selectedRow*radialDirection(1);
    Q(3*selectedDipoleIndex-1,:) = selectedRow*radialDirection(2);
    Q(3*selectedDipoleIndex,  :) = selectedRow*radialDirection(3);
    electrodePotentials = GainMat * Q;
end
figure;
timePoints = 1:size(electrodePotentials, 2);

for i = 1:numElctrodes
    subplot(7, 3, i);
    plot(electrodePotentials(i, :));
    title(ElectrodeNames{i});    
    xlim([timePoints(1), timePoints(end)]);
end

averagePotentials = zeros(1, numElctrodes);

for i = 1:numElctrodes 
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    averagePotentials(i) = electrode_peak_potential/(length(locations)*windowSize);
end

figure;
Display_Potential_3D(ModelParams.R(3), averagePotentials)

alpha = 1;
Q_MNE = GainMat'*inv(GainMat*GainMat'+eye(numElctrodes))*averagePotentials';

Omega = zeros(numDipoles, numDipoles);
for i = 1:numDipoles
    temp = 0;
    for j = 1:numElctrodes
        temp = temp + GainMat(j, 3*i-2:3*i)*GainMat(j, 3*i-2:3*i)';
    end
    Omega(i, i) = sqrt(temp);
end
W = kron(Omega, eye(3));
Wm = inv(W'*W);
alpha_w = 0.01;
Q_WMNE = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_w*eye(numElctrodes))*averagePotentials';

%% O(س)
clc;
Q_MNE_norms = zeros(numDipoles,1);
Q_WMNE_norms = zeros(numDipoles,1);

for i=1:numDipoles
    Q_MNE_norms(i) = norm(Q_MNE(3*i-2:3*i));
    Q_WMNE_norms(i) = norm(Q_WMNE(3*i-2:3*i));
end

%% P(ش)
clc;
labels = zeros(1, numDipoles);
for i = 1:length(selected_dipoles_indices)
    labels(selected_dipoles_indices(i)) = 1;
end
% Compute ROC
[FPR_MNE, TPR_MNE, ~, AUC_MNE] = perfcurve(labels, Q_MNE_norms, 1);
[FPR_wMNE, TPR_wMNE, ~, AUC_wMNE] = perfcurve(labels, Q_WMNE_norms, 1);

% Plot ROC curves
figure; hold on;
plot(FPR_MNE, TPR_MNE, 'b', 'LineWidth', 2);
plot(FPR_wMNE, TPR_wMNE, 'r', 'LineWidth', 2);
legend(['MNE (AUC = ' num2str(AUC_MNE) ')'], ['wMNE (AUC = ' num2str(AUC_wMNE) ')']);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curves for MNE and wMNE');
hold off;


%% Q(ص)
clc;
selected_dipoles_indices = [714, 715, 728, 729, 880, 881, 893, 894, 895, 896, 906, 907, 908, 909, 919, 920, 1050, 1051, 1062, 1063];
for i=1:length(selected_dipoles_indices)
    selectedRow = Interictal(i, :);
    selectedDipoleIndex = selected_dipoles_indices(i);
    radialDirection = LocMat(:, selectedDipoleIndex)/norm(LocMat(:, selectedDipoleIndex));
    Q = zeros(3*numDipoles, length(Interictal));
    Q(3*selectedDipoleIndex-2,:) = selectedRow*radialDirection(1);
    Q(3*selectedDipoleIndex-1,:) = selectedRow*radialDirection(2);
    Q(3*selectedDipoleIndex,  :) = selectedRow*radialDirection(3);
    electrodePotentials = GainMat * Q;
end

averagePotentials = zeros(1, numElctrodes);
for i = 1:numElctrodes 
    electrode_peak_potential = 0;
    for loc = locations
        elc_win = electrodePotentials(i, loc-halfWindow:loc+halfWindow);
        electrode_peak_potential = electrode_peak_potential + sum(elc_win);
    end
    averagePotentials(i) = electrode_peak_potential/(length(locations)*windowSize);
end

A1 = zeros(numDipoles, numDipoles);
minnorm = 1000;
for i=1:numDipoles
    for j=1:numDipoles
        if norm(LocMat(:,i)-LocMat(:,j)) == 1
            A1(i,j) = 1/6;
        end
    end
end
onep = ones(numDipoles,1);
A0 = inv(diag(A1*onep)+eye(numDipoles))*A1;
A = kron(A0, eye(3));
B = 6*(A-eye(3*numDipoles));
Omega = zeros(numDipoles, numDipoles);
for i = 1:numDipoles
    temp = 0;
    for j = 1:numElctrodes
        temp = temp + GainMat(j, 3*i-2:3*i)*GainMat(j, 3*i-2:3*i)';
    end
    Omega(i, i) = sqrt(temp);
end
Ww = kron(Omega, eye(3));
W = Ww * B' * B * Ww;
Wm = inv(W'*W);
alpha_loreta = 0.01;
Q_loreta_set = Wm*GainMat'*inv(GainMat*Wm*GainMat'+alpha_loreta*eye(numElctrodes))*averagePotentials';
CovarianceMatrix = electrodePotentials*electrodePotentials'/length(electrodePotentials);
Cin = inv(CovarianceMatrix);
alpha_sloreta = 0.01;
Q_sloreta_set = inv(GainMat'*Cin*GainMat+alpha_sloreta*eye(3*numDipoles))*GainMat'*Cin*averagePotentials';

%----------------------------------------
Q_loreta_norms = zeros(1,numDipoles);
Q_sloreta_norms = zeros(1,numDipoles);
for i=1:numDipoles
    Q_loreta_norms(i) = norm(Q_loreta_set(3*i-2:3*i));
    Q_sloreta_norms(i) = norm(Q_sloreta_set(3*i-2:3*i));
end

labels = zeros(1, numDipoles);
for i = 1:length(selected_dipoles_indices)
    labels(selected_dipoles_indices(i)) = 1;
end
% Compute ROC
[FPR_MNE, TPR_MNE, ~, AUC_MNE] = perfcurve(labels, Q_loreta_norms, 1);
[FPR_sloreta, TPR_sloreta, ~, AUC_sloreta] = perfcurve(labels, Q_sloreta_norms, 1);

% Plot ROC curves
figure; hold on;
plot(FPR_MNE, TPR_MNE, 'b', 'LineWidth', 2);
plot(FPR_sloreta, TPR_sloreta, 'r', 'LineWidth', 2);
legend(['Loreta (AUC = ' num2str(AUC_MNE) ')'], [['sLoreta (AUC = ' num2str(AUC_sloreta) ')']]);
xlabel('False Positive Rate');
ylabel('True Positive Rate');
title('ROC Curves for MNE and wMNE');
hold off;


