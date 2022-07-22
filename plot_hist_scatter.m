function plot_hist_scatter(performance, id, meth, save_flag)
% _
%   This function produces the scatterplot of responses and true directions.
%   The color represent the absolute deviation from the correct direction.
%
% IN: performance: a structure obtained by the function rdk_extract that
%     must include the cells directions,responses, and deviations with size 
%     length(s) X length(id.coh). These cells contain trialwise direction, 
%     responses and angular deviations (in degrees). The colormap of the
%     scatterplot indicate the absolute deviation of each datapoint.
%
%     id: is a structure obtained by the function rdk_extract that must
%     include the vector .st_id of size 1 x length(s) and contatining
%     the stimulus identity for every subject (1= transparent motion,
%     2=brownian motion, 3= white noise motion) and the eventual cell .rm
%     (size length(s) X length(id.coh)), containing the response map
%     information for every trial. If there is only one response map, then
%     all the entries in the cell will be = to the corresponding response
%     map number.
%
%     .coh is a vector with the coherence levels of the rdk shown during
%     experiment, sorted in descendent order.
%
%     s: represents the subject/s number/s. It can be an integer number or a
%     vector with multiple subjects.
%
%     rm: contains an integer number identifying the trials with a particular
%     response map. if rm is ~=0 then the analysis will be restricted to
%     the trials for which id.rm{k,j}==rm, otherwise all the trials will be
%     included in the analysis.
%
% OUT: one scatterplot for each stimulus identity with data pooled from all
%      the subjects in s
%
% written by Riccardo Barbieri <riccardo.barbieri@bccn-berlin.de>, 24/10/2018
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, ca. 2020
% finalized by Joram Soch <joram.soch@bccn-berlin.de>, 08/07/2022, 16:41
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 19/07/2022


% get stimuli and coherence level
% stim  = id.st_name;
stim    = {'TM', 'BM', 'WM'};
coh_str = {'100','50','25','12.5','0'};

% create scatter plots and histograms

% for each stimulus
for stim_ = 1:length(stim)

    % get subjects matching the current stimulus type
    ids   = id.st_id==stim_;
    cur_s = find(ids);

    % for each coherence level
    for coh_ = size(id.coh,2):-1:1

        all_dir{coh_}  = horzcat(performance.directions{cur_s,coh_});
        all_resp{coh_} = horzcat(performance.responses{cur_s,coh_});
        all_dev{coh_}  = horzcat(performance.deviations{cur_s,coh_});

        h = findobj('type','figure');
        n = length(h);
        figure_han = figure(n+1);
        figure_pos = [50+(coh_-1)*100, 50+(numel(stim)-stim_)*250-(coh_-1)*25, 620, 500];
        set(gcf, 'Position', figure_pos, 'Name', sprintf('Scatterplot of trial-wise responses (%s: %s%%)', stim{stim_}, coh_str{coh_}));
        
        create_scatterhist(all_dir{coh_}, all_resp{coh_}, all_dir{coh_}, all_resp{coh_}, figure_han, meth, stim{stim_}, coh_str{coh_}, save_flag)

    end
    
end