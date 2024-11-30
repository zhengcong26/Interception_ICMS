%% Script to create and initial 'unstable' weight matrix and a stability-optimised circuit

% This code creates an initial 'unstable' weight matrix W according to
% (Hennequin et al., Neuron, 2014) and then creates the stability-optimised
% variant so that the resulting neuronal dynamics display rich activity transients.
%
% Written by Jake Stroud
% modifiend by Cong Zheng

% Create intial 'unstable' weight matrix
N = 200;            %Number of neurons
p = 0.1;            %Density of connections
R = 10;             %Initial approximate spectral abscissa prior to stability optimisation
gamma = 3;          %The inhibition/excitation ratio
IEratio = 1/2;

W = initialnet_modified(N, p, R, gamma, IEratio); %Create initial 'unstable' weight matrix

% Create stabiliy-optimised circuit using the initial weight matrix W
rate = 10;              %Gradient-descent learning rate
desired_SA = 0.15;      %Ultimate desired spectral abscissa after stability optimisation

% Create stability-optimised circuit
Wsoc = soc_function(W, rate, desired_SA, gamma, IEratio);

%%
data=importdata('E:\Laplace task collection\Laplace code\Seismic_color_map.txt');
Seismic=data(:,1:3);

figure
set(gcf,'position',[30,300,300,230]);
h = heatmap(Wsoc);
% h = heatmap(W);
h.CellLabelColor = 'none';
h.GridVisible = 'off';
colormap(Seismic);
caxis([-4 4])
h.XDisplayLabels = nan(size(h.XDisplayData));
h.YDisplayLabels = nan(size(h.YDisplayData));
% axis square
% axis equal


