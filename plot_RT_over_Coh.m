function plot_RT_over_Coh(performance, id, meth, save_flag)
% _
% This function is written to plot the reaction times
% for the different stimuli and coherence levels.
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
% written by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 03/09/2020
% edited by Joram Soch <joram.soch@bccn-berlin.de>, 08/07/2022, 16:50
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 19/07/2022


color_code=['g','b','r'];

for stim_=1:length(id.st_name)
    
    if stim_==1
        figure('Name','Reaction time against coherence level','Color',[1 1 1],'Position',[50+0, 50+500, 600, 400])

        % JS
        hold on;
        plot(-1, -1, strcat('-',color_code(1)), 'LineWidth', 3);
        plot(-1, -1, strcat('-',color_code(2)), 'LineWidth', 3);
        plot(-1, -1, strcat('-',color_code(3)), 'LineWidth', 3);
        % JS
    end
    
    % CALCULATION
    
    rt=performance.rt(id.st_id==stim_,:);
    mean_rt= cellfun(@nanmean,rt);
    
    std_rt=std(mean_rt,[],1);
    [num,~]=size(mean_rt);
    se_rt=std_rt/sqrt(num);    
    
    mean_mean_rt = fliplr(mean(mean_rt,1));
    
    results(stim_).RT_COH=mean_mean_rt;
    results(stim_).RT_COH_SE=se_rt;
    
    % PLOTTING
    
    % plot with jitter
    if stim_==1
        delta =-0.08;
    elseif stim_==2
        delta = 0;
    elseif stim_==3
        delta = 0.08;
    end
    
    h=errorbar([1+delta 2+delta 3+delta 4+delta 5+delta],mean_mean_rt,se_rt,color_code(stim_),'LineStyle','none', 'LineWidth', 3);    
    h.CapSize = 8;
    % hold on  % JS
    plot([1+delta 2+delta 3+delta 4+delta 5+delta],mean_mean_rt,color_code(stim_) , 'LineWidth', 3)
    % hold off % JS
    hold on
    
end

% LAYOUT

    % JS
    set(gcf,'Position', [50, 50, 800, 450]);
    set(gcf,'Color', [1 1 1]);
    
    set(gca,'box','off');
    set(gca,'LineWidth',2);
    
    % x axis
    xlim([0.7 5.3])
    set(gca,'XTick',[1:1:5]);
    labels=char([{'0'};{'12.5'};{'25'};{'50'};{'100'}]);
    set(gca,'XTickLabel',labels, 'FontSize',30,'FontWeight','bold', 'LineWidth', [4]);
    xlabel('Coherence level [%]', 'Fontsize', 32)
    
    % y axis
    ylim([0.5,1.5]);    
    ylabel('Reaction time [s]', 'Fontsize', 32)
    legend('TM','BM','WM','Location','NorthEast')     
    set(gca,'YTick',[0.5:0.5:1.5]);
    a=get(gca);
       
    % save or not
    if save_flag==1
        save_folder=pwd;              
        save_name=[meth ' ' 'rt_over_coh'   '.pdf'];
        saveas(gcf,fullfile(save_folder,save_name))
    end
    
end