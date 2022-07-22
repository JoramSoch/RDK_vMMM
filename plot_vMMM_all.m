function plot_vMMM_all(performance, id, meth, save_flag)
% _
% Estimation of von Mises mixture model for random dot kinematogram data.
%
% written by Joram Soch <joram.soch@bccn-berlin.de>, 07/08/2018, 15:10 (V1)
% edited by Joram Soch <joram.soch@bccn-berlin.de>, 06/04/2020, 22:30 (V5)
% edited by Felix Töpfer <felix.toepfer@bccn-berlin.de>, 04/2020
% finalized by Joram Soch <joram.soch@bccn-berlin.de>, 08/07/2022, 17:45
% finalized by Joram Soch <joram.soch@bccn-berlin.de>, 22/07/2022, 11:38


% specify bias
bias   = 'model';               % bias is estimated as model parameter
% bias = 'preproc';             % bias is removed in pre-processing


%%% Step 0: data analysis preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
tasks = {'TM', 'BM', 'WM'};
clvls = id.coh;
level = {'100 %', '50 %', '25 %', '12.5 %', '0 %'};
model = {'m00', 'm01', 'm10', 'm11'};
color_code = ['g','b','r'];

% extract tasks
task_id = id.st_id;
dir_dev = performance.deviations;

% get numbers
num_task = max(task_id);
num_lvls = size(dir_dev,2);
num_subj = zeros(1,num_task);
num_trls = numel(dir_dev{1,1});
num_mods = numel(model);

% load data
raw_data = cell(max(num_subj),num_lvls,num_task);
for i = 1:num_task
    ind_subj    = find(task_id==i);
    num_subj(i) = numel(ind_subj);
    for j = 1:num_subj(i)
        for k = 1:num_lvls
            raw_data{j,k,i} = dir_dev{ind_subj(j),k};
        end;
    end;
end;
pp_data = cell(max(num_subj)+1,num_lvls,num_task);

% prepare results
m_bias = cell(max(num_subj)+1,num_lvls,num_task);
r_est  = cell(max(num_subj)+1,num_lvls,num_task);
m_est  = cell(max(num_subj)+1,num_lvls,num_task);
k_est  = cell(max(num_subj)+1,num_lvls,num_task);
r_avg  = cell(num_task);
m_avg  = cell(num_task);
k_avg  = cell(num_task);
r_SEs  = cell(num_task);
m_SEs  = cell(num_task);
k_SEs  = cell(num_task);
MLL    = cell(max(num_subj)+1,num_lvls,num_task);
AIC    = cell(max(num_subj)+1,num_lvls,num_task);
AICw   = cell(max(num_subj)+1,num_lvls,num_task);


% Task 1: model estimation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fprintf('\n');

% for every task/condition/stimulus
for i = 1:num_task
    
    fprintf('-> Stimulus %d: ', i);
    
    % preallocate averages
    r_avg{i} = zeros(num_lvls, 3);                      % avg over subjects
    m_avg{i} = zeros(num_lvls, 2);
    k_avg{i} = zeros(num_lvls, 2);
    r_SEs{i} = zeros(num_lvls, 3);                      % SE over subjects
    m_SEs{i} = zeros(num_lvls, 2);
    k_SEs{i} = zeros(num_lvls, 2);
    
    % for every subject
    for j = 1:num_subj(i)+1
        
        if j <= num_subj(i), fprintf('%d, ', j); end;
        if j >  num_subj(i), fprintf('all, ');   end;
        
        % pre-process each subject
        if j <= num_subj(i)
            % remove bias in pre-processing
            if strcmp(bias,'preproc')
                % JS, vMMM V4
                % [pp_data{i}{j,1}, m_bias{i}(j,1)] = ME_preproc(raw_data{i}{j,1});
                % for k = 2:num_lvls
                %     pp_data{i}{j,k} = ME_preproc(raw_data{i}{j,k}, m_bias{i}(j,1));
                % end;
                
                % JS, vMMM V5
                % [pp_data{j,1,i}, m_bias{j,1,i}] = ME_preproc(raw_data{j,1,i});
                % for k = 2:num_lvls
                %     pp_data{j,k,i} = ME_preproc(raw_data{j,k,i}, m_bias{j,1,i});
                % end;
                
                % FMT, 08/04/2020
                m_bias{j,1,i} = calc_bias(squeeze(raw_data(j,:,i)));
                for k = 1:num_lvls
                    pp_data{j,k,i} = ME_preproc(raw_data{j,k,i}, m_bias{j,1,i});
                end;
            end;
            % do not remove the bias
            if strcmp(bias,'model')
                for k = 1:num_lvls
                    [pp_data{j,k,i}, m_bias{j,k,i}] = ME_preproc(raw_data{j,k,i}, [], false);
                end;
            end;
        end;
        
        % concatenate all subjects
        if j > num_subj(i)
            for k = 1:num_lvls
                y = [];
                for jc = 1:num_subj(i)
                    y = [y; pp_data{jc,k,i}];
                end;
                pp_data{j,k,i} = y;
            end;
        end;
        
        % analyze subject's data
        for k = 1:num_lvls
            r_est{j,k,i} = zeros(num_mods+1, 3);
            m_est{j,k,i} = zeros(num_mods+1, 2);
            k_est{j,k,i} = zeros(num_mods+1, 2);
            MLL{j,k,i}   = zeros(num_mods, 1);
            for l = 1:num_mods
                % exclude bias from modelling
                if strcmp(bias,'preproc')
                    [r_est{j,k,i}(l,:), m_est{j,k,i}(l,:), k_est{j,k,i}(l,:), MLL{j,k,i}(l)] = ME_vMMM_ML(pp_data{j,k,i}, [], [], model{l});
                end;
                % account for bias in model
                if strcmp(bias,'model')
                    [r_est{j,k,i}(l,:), m_est{j,k,i}(l,:), k_est{j,k,i}(l,:), MLL{j,k,i}(l)] = ME_vMMM_ML(pp_data{j,k,i}, [], [], strcat(model{l},'b'));
                end;
            end;
            % perform model averaging
            if strcmp(bias,'preproc'), p = [0, 2, 2, 4]'; end;
            if strcmp(bias,'model'),   p = [0, 3, 3, 5]'; end;
            AIC{j,k,i}  = -2*MLL{j,k,i} + 2*p;
            AICw{j,k,i} = exp(-1/2*(AIC{j,k,i} - min(AIC{j,k,i}))) ./ sum( exp(-1/2*(AIC{j,k,i} - min(AIC{j,k,i}))) );
            r_est{j,k,i}(num_mods+1,:) = sum( r_est{j,k,i}(1:num_mods,:) .* repmat(AICw{j,k,i},[1 size(r_est{j,k,i},2)]) );
            m_est{j,k,i}(num_mods+1,:) = sum( m_est{j,k,i}(1:num_mods,:) .* repmat(AICw{j,k,i},[1 size(m_est{j,k,i},2)]) );
            k_est{j,k,i}(num_mods+1,:) = sum( k_est{j,k,i}(1:num_mods,:) .* repmat(AICw{j,k,i},[1 size(k_est{j,k,i},2)]) );
        end;
        
    end;
    
    fprintf('done.\n');
    
    % for every coherence level
    r_hat = zeros(num_lvls, 3, num_subj(i)+1);
    m_hat = zeros(num_lvls, 2, num_subj(i)+1);
    k_hat = zeros(num_lvls, 2, num_subj(i)+1);
    for j = 1:num_subj(i)+1
        for k = 1:num_lvls
            r_hat(k,:,j) = r_est{j,k,i}(num_mods+1,:);
            m_hat(k,:,j) = m_est{j,k,i}(num_mods+1,:);
            k_hat(k,:,j) = k_est{j,k,i}(num_mods+1,:);
        end;
    end;
    r_avg{i} = mean(r_hat(:,:,1:end-1),3);
    m_avg{i} = mean(m_hat(:,:,1:end-1),3);
    k_avg{i} = mean(k_hat(:,:,1:end-1),3);
    r_SEs{i} = sqrt(var(r_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    m_SEs{i} = sqrt(var(m_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    k_SEs{i} = sqrt(var(k_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    
end;

fprintf('\n');


%%% Task 2: model comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%% Figure 9/10: Model Comparison (JS, 07/2022)

figure('Name','vMMM: model comparison','Color',[1 1 1],'Position',[50 50 1600 900]);
for i = 1:num_task
    
    PPs = zeros(num_lvls, num_mods);
    
    for k = 1:num_lvls
        
        PPs(k,:) = AICw{num_subj(i)+1,k,i}';
        
        subplot(num_task,num_lvls,(i-1)*num_lvls+(num_lvls-k+1));
        hold on;
        bar(PPs(k,:), 0.5, 'FaceColor', [0.5 0.5 0.5]);
        title(sprintf('%s: %s%', tasks{i}, level{k}));
        set(gca, 'FontSize', 15, 'FontWeight', 'bold', 'LineWidth', 4);
        xlim([-0.2 5.5])
        xticks([1:num_mods])
        xticklabels(model)
        xtickangle(315)
        xlabel('Models', 'FontSize', 15);

        ylim([0 1.1])
        yticks([0 0.5 1])
        yticklabels({'0','0.5','1'})
        ylabel('AIC weights', 'FontSize', 15);
        
    end

end


%%% Task 3: model parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Do you want to plot legends?
plot_legend = 1;
    
% specify plotting
axislabelsize = 32;
axisticksize  = 30;
axislinewidth = 4;
plotlinewidth = 3;


%%% Figure 11: Detection Frequency (FMT, 10/2019)

figure('Name','vMMM: detection frequency','Color',[1 1 1],'Position',[50+0, 50+500, 600, 400])
for i = 1:num_task
    r_avg_stim=r_avg{i}';
    r_se_stim=r_SEs{i}';
    
    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta 5+delta],fliplr(r_avg_stim(1,:)),fliplr(r_se_stim(1,:)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on
       
     results(i).detect_freq=fliplr(r_avg_stim(1,:));
     results(i).detect_freq_se=fliplr(r_se_stim(1,:));
     
end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
        % x axis
        xlim([0.7 5.3])
        set(gca,'XTick',[1:1:5]);
        labels=char([{'0'};{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        % y axis
        ylabel('Detection frequency', 'FontSize', axislabelsize)
        ylim([0,1]);

        if plot_legend
            legend('TM','BM','WM','Location','southeast')
        end
    
        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'detection frequency'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end         
        

%%% Figure 14: vMMM: ROOD frequency (FMT, 10/2019)

figure('Name','vMMM: ROOD frequency','Color',[1 1 1],'Position',[50+0, 50+250, 600, 400])
for i = 1:num_task
    r_avg_stim=r_avg{i}';
    r_se_stim=r_SEs{i}';

    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta 5+delta],fliplr(r_avg_stim(2,:)),fliplr(r_se_stim(2,:)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on
     
     results(i).rood_freq=fliplr(r_avg_stim(2,:));
     results(i).rood_freq_se=fliplr(r_se_stim(2,:));

end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
        % x axis
        xlim([0.7 5.3])
        set(gca,'XTick',[1:1:5]);
        labels=char([{'0'};{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        % y axis
        ylabel('ROOD frequency', 'FontSize', axislabelsize)
        ylim([0,0.4]);

        if plot_legend
            legend('TM','BM','WM','Location','northeast')
        end

        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'ROOD frequency'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end                 


%%% Figure 16: Guessing Frequency (FMT, 10/2019)

figure('Name','vMMM: guessing frequency','Color',[1 1 1],'Position',[50+0, 50+0, 600, 400])
for i = 1:num_task
    r_avg_stim=r_avg{i}';
    r_se_stim=r_SEs{i}';

       
    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta 5+delta],fliplr(r_avg_stim(3,:)),fliplr(r_se_stim(3,:)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on

     results(i).guess_rate=fliplr(r_avg_stim(3,:));
     results(i).guess_rate_se=fliplr(r_se_stim(3,:));

 
end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
        % x axis
        xlim([0.7 5.3])
        set(gca,'XTick',[1:1:5]);
        labels=char([{'0'};{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        % y axis
        ylabel('Guess frequency', 'FontSize',  axislabelsize)
        ylim([0,1]);

        if plot_legend
            legend('TM','BM','WM','Location','northeast')
        end
        
        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'guess frequency'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end                 
        

%%% Figure 12: Detection Precision (FMT, 10/2019)

figure('Name','vMMM: detection precision','Color',[1 1 1],'Position',[50+200, 50+500-25, 600, 400])
for i = 1:num_task
    k_avg_stim=k_avg{i}';
    k_se_stim=k_SEs{i}';
       
    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta],fliplr(k_avg_stim(1,1:4)),fliplr(k_se_stim(1,1:4)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on     
     
     results(i).detect_precision=fliplr(k_avg_stim(1,:));
     results(i).detect_precision_se=fliplr(k_se_stim(1,:));

end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
        % x axis
        xlim([0.7 4.3])
        set(gca,'XTick',[1:1:4]);
        labels=char([{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        % y axis
        ylabel('Detection precision', 'FontSize', axislabelsize)
        ylim([0,35]);
        
        if plot_legend
            legend('TM','BM','WM','Location','southeast')
        end
         
        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'detection precision'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end                 
        
        
%%% Figure 15: ROOD Precision (FMT, 10/2019)

figure('Name','vMMM: ROOD precision','Color',[1 1 1],'Position',[50+200, 50+250-25, 600, 400])
for i = 1:num_task
    k_avg_stim=k_avg{i}';
    k_se_stim=k_SEs{i}';
       
    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta],fliplr(k_avg_stim(2,1:4)),fliplr(k_se_stim(2,1:4)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on          

     results(i).rood_precision=fliplr(k_avg_stim(2,:));
     results(i).rood_precision_se=fliplr(k_se_stim(2,:));

end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
            % x axis
        xlim([0.7 4.3])
        set(gca,'XTick',[1:1:4]);
        labels=char([{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        % y axis
        ylabel('ROOD precision', 'FontSize',  axislabelsize)
        ylim([0,15]);

        if plot_legend
            legend('TM','BM','WM','Location','north')
        end

        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'ROOD precision'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end                 
        
        
%%% Figure 13: Detection Location (FMT, 10/2019)

figure('Name','vMMM: detection bias','Color',[1 1 1],'Position',[50+400, 50+500-50, 600, 400])
for i = 1:num_task
    m_avg_stim=rad2deg(m_avg{i}');
    m_se_stim=rad2deg(m_SEs{i}');

    if i==1
        delta =-0.08;
    elseif i==2
        delta = 0;
    elseif i==3
        delta = 0.08;
    end
    
        h=errorbar([1+delta 2+delta 3+delta 4+delta],fliplr(m_avg_stim(1,1:4)),fliplr(m_se_stim(1,1:4)),color_code(i), 'LineWidth', plotlinewidth);    
        h.CapSize = 8;
        hold on          
     
     results(i).detect_precision=fliplr(m_avg_stim(1,:));
     results(i).detect_precision_se=fliplr(m_se_stim(1,:));

end    
            
        set(gcf,'color','w');
        set(gca,'box','off');
        set(gca,'LineWidth',axislinewidth);
        
            % x axis
        xlim([0.7 4.3])
        set(gca,'XTick',[1:1:4]);
        labels=char([{'12.5'};{'25'};{'50'};{'100'}]);
        set(gca,'XTickLabel',labels, 'FontSize',axisticksize,'FontWeight','bold', 'LineWidth', axislinewidth);
        xlabel('Coherence level [%]', 'Fontsize', axislabelsize);
        
        
        % y axis
        ylabel('Detection mean [°]', 'FontSize', axislabelsize)
        ylim([-10,20]);
        yline(0,'--', 'LineWidth',plotlinewidth);

        if plot_legend
            legend('TM','BM','WM','Location','southeast')
        end
              
        if save_flag==1
            save_folder=pwd;              
            save_name=[meth ' ' 'detection location'   '.pdf'];
            saveas(gcf,fullfile(save_folder,save_name))
        end