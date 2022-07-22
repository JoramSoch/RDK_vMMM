% RDK-vMMM: behavioral analysis
% _
% This script reproduces the analysis pipeline from Töpfer et al.
% (in review at JoV) and creates figures presented in the paper.
% 
% written by Felix Töpfer <felix.toepfer@bccn-berlin.de>, ca. 2020
% finalized by Joram Soch <joram.soch@bccn-berlin.de>, 08/07/2022, 16:01

clear
close all


%%% load data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% add paths
curr_dir = pwd;
addpath(strcat(curr_dir,'/tools/'));
addpath(strcat(curr_dir,'/tools/vMMM/'));

% Which response method you want to analyze? 
meth = 'trackball';
% 'bar' for rotating bar or
% 'trackball' for trackball

% rotating bar
if strcmp(meth,'bar')
    load(fullfile(curr_dir,'data','data_set_bar.mat'));
    performance = data_set_bar.performance;
    id = data_set_bar.id;

% trackball
elseif strcmp(meth,'trackball')
    load(fullfile(curr_dir,'data','data_set_trackball.mat'));
    performance = data_set_trackball.performance;
    id = data_set_trackball.id;

% error message
else
    error('No such response method. Please enter ''bar'' or ''trackball''.')
    
end;


%%% create figures  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Figure 3/4: stimulus and response distributions
plot_hist_scatter(performance, id, meth, 0);

% Table 1/2: Cramér-von-Mises test statistics
stats = calc_Cramer_vM(performance, id);

% Figure 5A/B: reaction time over coherence
plot_RT_over_Coh(performance, id, meth, 0);

% Figure 6/7: histogram of response deviations
plot_hist_deviation(performance, id, meth, 0)

% Figures 9-16: estiamtes from vMMM
plot_vMMM_all(performance, id, meth, 0);