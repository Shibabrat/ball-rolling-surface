axesFontName = 'factory';
%     axesFontName = 'Times New Roman';
set(0,'Defaulttextinterpreter','latex', ...
    'DefaultAxesFontName', axesFontName, ...
    'DefaultTextFontName', axesFontName, ...
    'DefaultAxesFontSize',16, ...
    'DefaultTextFontSize',16, ...
    'Defaultuicontrolfontweight','bold', ...
    'Defaulttextfontweight','bold', ...
    'Defaultaxesfontweight','bold');

load mfdU1p_stable_L12_energy_case1.mat %<====check if data file is in the path

energyTubeTraj = get_energy_points_ball_rolling(xW);
plot(energyTubeTraj-energyTube)

xlabel('$\#$ trajectories','fontsize',18);
ylabel('$E_{trajectories} - E_{target}$','fontsize',18)

