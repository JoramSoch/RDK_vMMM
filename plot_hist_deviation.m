function plot_hist_deviation(performance, id, meth, save_flag)
% _
% This function is written to plot the deviation (theta = physical stimulus 
% direction - recorded response) of responses from the presented motion
% direction for every stimulus and coherence pooled over all subjects
% 
% IN:   performance:    structre, containing the recorded behavior for every
%                       subject
%       id:             structure, identification structure that links the
%                       recorded behavior in the structure performance to
%                       the coherence level, stimulus type,...
%       meth:           string, the name of the response method you want to
%                       plot
%       save_flag:      number, can take 0 (= do not save the created 
%                       figures) and 1 (save created figures)
% 
% written by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 02/09/2020
% finalized by Joram Soch <joram.soch@bccn-berlin.de>, 08/07/2022, 16:55
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 11/06/2022
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 19/07/2022


% for each stimulus

figure('Name','Histogram of trial-wise response deviations','Color',[1 1 1],'Position',[50 50 1600 900]);

for stim_ = 1:length(id.st_name)

    % get subjects matching the current stimulus type
    ids   = id.st_id==stim_;
    cur_s = find(ids);
    tasks = {'TM', 'BM', 'WM'};
    level = {'100%', '50%', '25%', '12.5%', '0%'};

    % for each coherence level
    for coh_ = 1:size(id.coh,2)

        dev = [];        
        dev = horzcat(performance.deviations{cur_s,coh_});
        dev = dev(~isnan(dev));
        dev = deg2rad(dev);

        subplot(length(id.st_name),size(id.coh,2),(stim_-1)*size(id.coh,2)+(size(id.coh,2)-coh_+1));
        hold on;
        
        histogram(dev,60,'FaceColor',[0 0 0])
        axis([-pi, +pi, 0, 250]);
        xline(0, '--', 'LineWidth', 2);
        set(gca,'Box','On');
        set(gca, 'XTick', [-pi:(pi/2):+pi], 'XTickLabel', [-180:90:+180], 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', [2]);
        xlabel('\theta [°]', 'FontSize', 15);
        ylabel('counts', 'FontSize', 15);
        yticks([0 125 250])        
        title_str = sprintf('%s: %s%', tasks{stim_}, level{coh_});
        title(title_str, 'FontSize', 18, 'FontWeight', 'bold');

        if save_flag == 1
            save_folder = pwd;              
            save_name   = [meth ' ' 'hist_deviation'  ' ' id.st_name{stim_} ' ' mat2str(id.coh(1,coh_)*100) '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end

    end
    
end