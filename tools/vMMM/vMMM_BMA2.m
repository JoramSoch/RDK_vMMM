% von Mises Mixture Model for Random Dot Motion
% _
% This script performs a vMMM analysis for data from the RDM experiment.
%
% Author: Joram Soch, BCCN Berlin
% E-Mail: joram.soch@bccn-berlin.de
% 
% First edit: 07/08/2018, 15:10 (V1)
%  Last edit: 07/04/2020, 16:50 (V5)


clear
close all

% specify bias
bias   = 'model';               % bias is estimated as model parameter
% bias = 'preproc';             % bias is removed in pre-processing


%%% Step 0: data analysis preparation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data
load data_36_subj.mat
tasks = id.t_name;
clvls = id.coh;
level = {'100 %', '50 %', '25 %', '12.5 %', '0 %'};
model = {'m00', 'm01', 'm10', 'm11'};

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
raw_data = cell(max(num_subj),num_task);
for i = 1:num_task
    ind_subj    = find(task_id==i);
    num_subj(i) = numel(ind_subj);
    for j = 1:num_subj(i)
        raw_data{j,i} = cell(num_lvls,1);
        for k = 1:num_lvls
            raw_data{j,i}{k} = dir_dev{ind_subj(j),k};
        end;
    end;
end;
pp_data = cell(max(num_subj)+1,num_task);

% prepare results
m_bias = cell(max(num_subj)+1,num_task);
r_est  = cell(max(num_subj)+1,num_task);
m_est  = cell(max(num_subj)+1,num_task);
k_est  = cell(max(num_subj)+1,num_task);
r_avg  = cell(num_task);
m_avg  = cell(num_task);
k_avg  = cell(num_task);
r_SEs  = cell(num_task);
m_SEs  = cell(num_task);
k_SEs  = cell(num_task);
AICw   = cell(max(num_subj)+1,num_task);


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
        
        % preallocate pre-processing
        pp_data{j,i} = cell(num_lvls, 1);
        m_bias{j,i}  = zeros(num_lvls, 1);
        
        % pre-process each subject
        if j <= num_subj(i)
            % remove bias in pre-processing
            if strcmp(bias,'preproc')
                [pp_data{j,i}{1}, m_bias{j,i}(1)] = ME_preproc(raw_data{j,i}{1});
                for k = 2:num_lvls
                    pp_data{j,i}{k} = ME_preproc(raw_data{j,i}{k}, m_bias{j,i}(1));
                end;
            end;
            % do not remove the bias
            if strcmp(bias,'model')
                for k = 1:num_lvls
                    [pp_data{j,i}{k}, m_bias{j,i}(k)] = ME_preproc(raw_data{j,i}{k}, [], false);
                end;
            end;
        end;
        
        % concatenate all subjects
        if j > num_subj(i)
            for k = 1:num_lvls
                y = [];
                for jc = 1:num_subj(i)
                    y = [y; pp_data{jc,i}{k}];
                end;
                pp_data{j,i}{k} = y;
            end;
        end;
        
        % analyze subject's data
        r_est{j,i} = zeros(num_lvls, 3);
        m_est{j,i} = zeros(num_lvls, 2);
        k_est{j,i} = zeros(num_lvls, 2);
        AICw{j,i}  = zeros(num_mods, num_lvls);
        for k = 1:num_lvls
            % exclude bias from modelling
            if strcmp(bias,'preproc')
                s = {'m00', 'm01', 'm10', 'm11'};
            end;
            % account for bias in model
            if strcmp(bias,'model')
                s = {'m00b', 'm01b', 'm10b', 'm11b'};
            end;
            [r_est{j,i}(k,:), m_est{j,i}(k,:), k_est{j,i}(k,:), AICw{j,i}(:,k)] = ME_vMMM_MA(pp_data{j,i}{k}, [], [], s);
        end;
        
    end;
    
    fprintf('done.\n');
    
    % for every coherence level
    r_hat = zeros(num_lvls, 3, num_subj(i)+1);
    m_hat = zeros(num_lvls, 2, num_subj(i)+1);
    k_hat = zeros(num_lvls, 2, num_subj(i)+1);
    for j = 1:num_subj(i)+1
        r_hat(:,:,j) = r_est{j,i};
        m_hat(:,:,j) = m_est{j,i};
        k_hat(:,:,j) = k_est{j,i};
    end;
    r_avg{i} = mean(r_hat(:,:,1:end-1),3);
    m_avg{i} = mean(m_hat(:,:,1:end-1),3);
    k_avg{i} = mean(k_hat(:,:,1:end-1),3);
    r_SEs{i} = sqrt(var(r_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    m_SEs{i} = sqrt(var(m_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    k_SEs{i} = sqrt(var(k_hat(:,:,1:end-1),[],3))./sqrt(num_subj(i));
    
end;

fprintf('\n');

% open figure (1): pooling, all
figure('Name','Model estimation: pooling (all)','Color',[1 1 1],'Position',[50 50 1600 900]);

% plot results
for i = 1:num_task
    
    % intervals
    dl  = ([1:num_lvls]-(num_lvls+1)/2)./(num_lvls+1.5);
    ddl = 1/(4*(num_lvls+1.5));
    
	% estimates
    r_hat = zeros(num_lvls, 3);
    m_hat = zeros(num_lvls, 2);
    k_hat = zeros(num_lvls, 2);
    for k = 1:num_lvls
        r_hat(k,:) = r_est{num_subj(i)+1,i}(k,:);
        m_hat(k,:) = m_est{num_subj(i)+1,i}(k,:);
        k_hat(k,:) = k_est{num_subj(i)+1,i}(k,:);
    end;
    
    % frequencies
    subplot(num_task,6,(i-1)*6+[1:3]);
    hold on;
    bar(r_hat', 'grouped');
    axis([(1-0.5), (3+0.5), 0, 1]);
    set(gca,'Box','On');
    if i == 3, legend(level, 'Location', 'North'); end;
    set(gca,'XTick',[1:3],'XTickLabel',{'r1 [detection]', 'r2 [ROOD]', 'r3 [guessing]'});
    ylabel(sprintf('%s (N = %d)', tasks{i}, num_subj(i)), 'FontSize', 20);
    if i == 1, title('Frequencies', 'FontSize', 20); end;
    
    % means
    subplot(num_task,6,(i-1)*6+[4]);
    hold on;
    bar(m_hat', 'grouped');
    axis([(1-0.5), (1+0.5), -(pi/6), +(pi/6)]);
    set(gca,'Box','On');
    set(gca,'XTick',[1],'XTickLabel',{'m1 [detection]'});
    set(gca,'YTick',[-(pi/6):(pi/18):+(pi/6)],'YTickLabel',[-30:10:+30]);
    if i == 1, title('Locations', 'FontSize', 20); end;
    
    % precision
    subplot(num_task,6,(i-1)*6+[5:6]);
    hold on;
    bar(k_hat', 'grouped');
    axis([(1-0.5), (2+0.5), 0, (15/10)*max(max(k_hat))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:2],'XTickLabel',{'k1 [detection]', 'k2 [ROOD]'});
    if i == 1, title('Precisions', 'FontSize', 20); end;
    
end;

% open figure (2): averaging, all
figure('Name','Model estimation: averaging (all)','Color',[1 1 1],'Position',[50 50 1600 900]);

% plot results
for i = 1:num_task
    
    % intervals
    dl  = ([1:num_lvls]-(num_lvls+1)/2)./(num_lvls+1.5);
    ddl = 1/(4*(num_lvls+1.5));
    
    % frequencies
    subplot(num_task,6,(i-1)*6+[1:3]);
    hold on;
    bar(r_avg{i}', 'grouped');
    for k = 1:num_lvls % SEs
        for l = 1:3
            fig_draw_CI(l+dl(k), r_avg{i}(k,l)+r_SEs{i}(k,l), r_avg{i}(k,l)-r_SEs{i}(k,l), 2*ddl, '-k', 2);
        end;
    end;
    axis([(1-0.5), (3+0.5), 0, 1]);
    set(gca,'Box','On');
    if i == 3, legend(level, 'Location', 'North'); end;
    set(gca,'XTick',[1:3],'XTickLabel',{'r1 [detection]', 'r2 [ROOD]', 'r3 [guessing]'});
    ylabel(sprintf('%s (N = %d)', tasks{i}, num_subj(i)), 'FontSize', 20);
    if i == 1, title('Frequencies', 'FontSize', 20); end;
    
    % means
    subplot(num_task,6,(i-1)*6+[4]);
    hold on;
    bar(m_avg{i}', 'grouped');
    for k = 1:num_lvls % SEs
        for l = 1:1
            fig_draw_CI(l+dl(k), m_avg{i}(k,l)+m_SEs{i}(k,l), m_avg{i}(k,l)-m_SEs{i}(k,l), 2*ddl, '-k', 2);
        end;
    end;
    axis([(1-0.5), (1+0.5), -(pi/6), +(pi/6)]);
    set(gca,'Box','On');
    set(gca,'XTick',[1],'XTickLabel',{'m1 [detection]'});
    set(gca,'YTick',[-(pi/6):(pi/18):+(pi/6)],'YTickLabel',[-30:10:+30]);
    if i == 1, title('Locations', 'FontSize', 20); end;
    
    % precision
    subplot(num_task,6,(i-1)*6+[5:6]);
    hold on;
    bar(k_avg{i}', 'grouped');
    for k = 1:num_lvls % SEs
        for l = 1:2
            fig_draw_CI(l+dl(k), k_avg{i}(k,l)+k_SEs{i}(k,l), k_avg{i}(k,l)-k_SEs{i}(k,l), 2*ddl, '-k', 2);
        end;
    end;
    axis([(1-0.5), (2+0.5), 0, (15/10)*max(max(k_avg{i}))]);
    set(gca,'Box','On');
    set(gca,'XTick',[1:2],'XTickLabel',{'k1 [detection]', 'k2 [ROOD]'});
    if i == 1, title('Precisions', 'FontSize', 20); end;
    
end;


%%% Task 2: model comparison %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% open figure (3): pooling
figure('Name','Model comparison: pooling','Color',[1 1 1],'Position',[50 50 900 900]);

% plot results
for i = 1:num_task
    % Akaike weights
    subplot(num_task,1,i);
    bar(AICw{num_subj(i)+1,i}, 'grouped');
    axis([(1-0.5), (4+0.5), 0, 1]);
    set(gca,'Box','On');
    if i == 1, legend(level, 'Location', 'NorthWest'); end;
    set(gca,'XTick',[1:4],'XTickLabel',{'m00 [only guessing]', 'm01 [guessing + ROOD]', 'm10 [guessing + detection]', 'm11 [guessing + ROOD + detection]'});
    ylabel(sprintf('%s (N = %d)', tasks{i}, num_subj(i)), 'FontSize', 20);
    if i == 1, title('AIC weights', 'FontSize', 20); end;
end;

% open figure (4): averaging
figure('Name','Model comparison: separate','Color',[1 1 1],'Position',[50 50 1800 900]);

% plot results
for i = 1:num_task
    for j = 1:num_subj(i)
        % Akaike weights
        subplot(num_task,max(num_subj),(i-1)*max(num_subj)+j);
        bar(AICw{j,i}, 'grouped');
        axis([(1-0.5), (4+0.5), 0, 1]);
        set(gca,'Box','On');
        set(gca,'XTick',[1:4],'XTickLabel',{'m00', 'm01', 'm10', 'm11'});
        if j == 1
            ylabel(sprintf('%s (N = %d)', tasks{i}, num_subj(i)), 'FontSize', 20);
            if i == 1, title('AIC weights', 'FontSize', 20); end;
        end;
    end;
end;