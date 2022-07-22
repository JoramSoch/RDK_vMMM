function [stats]=calc_Cramer_vM(performance,id)



stim={'Trans','Brown','White'};

coh_str=({'100','50','25','12.5','0'});


% allocate

stats.cramer_vM=cell(3,5);

%

for stim_=1:length(stim)


    ids=id.st_id==stim_;%get subjects matching the current stimulus type
    cur_s=find(ids);

    for coh_=1:size(id.coh,2)

        all_dir{coh_}=horzcat(performance.directions{cur_s,coh_});
        all_resp{coh_}=horzcat(performance.responses{cur_s,coh_});
        all_dev{coh_}=horzcat(performance.deviations{cur_s,coh_});

%         h =  findobj('type','figure');
%         n = length(h);
%         figure_handel = figure(n+1);
% 
%         create_scatterhist(all_dir{coh_},all_resp{coh_},all_dir{coh_},all_resp{coh_},figure_handel, meth, stim{stim_}, coh_str{coh_}, save_flag)

        
        [h0,p,cmstats]= cmtest2(all_dir{coh_},all_resp{coh_});
        stats.cramer_vM{stim_,coh_} = [h0,p,cmstats];
        
    end
    % flip left to right that it has the same order as the figures 

end

    stats.cramer_vM=fliplr(stats.cramer_vM);





end

