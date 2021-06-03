%% ============= MAIN FUNCTION ============== %%
% ----------------------------------------------
% Parameters Set Up 
% ----------------------------------------------
COD_file_path        = 'COD_sample.cod';
start_time           = 0;    % in sec
end_time             = 5400; % in sec
scan_duration_in_use = end_time - start_time;
alpha                = 1.4;  % Options: 1,1.2,1.3,1.4,1.6,1.8,2,etc.
n_max                = 300;  % Recommended for a 90-min scan
Scout_segments       = 6;    % The number in which the scan will be divided into

% ----------------------------------------------
% Calculating motion time points (MTPs)
% ----------------------------------------------
cod_direction_pool        = {'x','y','z'};
for idx_cod_dir           = 1:length(cod_direction_pool)
    cod_direction         = cod_direction_pool{idx_cod_dir};
    MTP                   = Get_MTPs(COD_file_path,start_time,end_time,alpha,n_max,Scout_segments,cod_direction);
    MTP_pool{idx_cod_dir} = MTP;
end
MTP_x = MTP_pool{1};
MTP_y = MTP_pool{2};
MTP_z = MTP_pool{3};

% ----------------------------------------------
% Combine the MTPs from all 3 COD directions
% ----------------------------------------------
union1       = union(MTP_x,MTP_y);
union_all    = union(union1,MTP_z);
MFF_start    = [1;union_all];
MFF_end      = [union_all-1;scan_duration_in_use];
MFF_duration = MFF_end-MFF_start+1;

for idx_MFF               = 1:numel(union_all)+1
    MFF(idx_MFF).start    = MFF_start(idx_MFF);
    MFF(idx_MFF).end      = MFF_end(idx_MFF);
    MFF(idx_MFF).duration = MFF_duration(idx_MFF);
end
save('MFF_union_all.mat','MFF');

% ----------------------------------------------
% Plotting COD traces with MTPs
% ----------------------------------------------
figure; 
for idx_cod_dir   = 1:length(cod_direction_pool)
    cod_direction = cod_direction_pool{idx_cod_dir};
    allcod        = dlmread(COD_file_path,' ',1,0);
    
    if strcmp(cod_direction,'x'); cod = allcod(:,2); end
    if strcmp(cod_direction,'y'); cod = allcod(:,3); end
    if strcmp(cod_direction,'z'); cod = allcod(:,4); end
    
    cod           = cod(:,1);
    prompt        = allcod(:,end);
    cod_in_use    = cod(start_time+1:end_time,1);
    prompt_in_use = prompt(start_time+1:end_time,1);
    
    for idx_MFF = 1:numel(union_all)+1
        MFF_each_dir(idx_MFF).start       = MFF_start(idx_MFF);
        MFF_each_dir(idx_MFF).end         = MFF_end(idx_MFF);
        MFF_each_dir(idx_MFF).duration    = MFF_duration(idx_MFF);
        MFF_each_dir(idx_MFF).mean        = mean(cod_in_use(MFF_each_dir(idx_MFF).start:MFF_each_dir(idx_MFF).end));
        MFF_each_dir(idx_MFF).std         = std(cod_in_use(MFF_each_dir(idx_MFF).start:MFF_each_dir(idx_MFF).end));
        MFF_each_dir(idx_MFF).residue     = MFF_duration(idx_MFF)*MFF_each_dir(idx_MFF).std.^2;
        MFF_each_dir(idx_MFF).prompt_rate = sum(prompt_in_use(MFF_each_dir(idx_MFF).start:MFF_each_dir(idx_MFF).end))/(MFF_each_dir(idx_MFF).end-MFF_each_dir(idx_MFF).start+1);
    end
    save(['COD ',cod_direction,'_MFF_union_all.mat'],'MFF_each_dir');
    
    % ----- For Visualization: overlaying MTPs on the COD trace to show MFFs
    MFF_indication = zeros(numel(cod_in_use),1)+min(cod_in_use);
    for idx = 1:numel(union_all)
        MFF_indication(MFF_each_dir(idx).start:MFF_each_dir(idx).end-1)=1*max(cod_in_use);
    end
    
    % ----- Plotting MTPs in 3 COD directions
    subplot(3,1,idx_cod_dir)
    plot(cod_in_use);hold on;
    plot(MFF_indication);hold on;
    title({['COD',cod_direction,' - MFF\_num = ',int2str(numel(union_all)),...
        ',inflation\_factor = ',num2str(alpha),...
        ',start\_time = ',int2str(start_time)]});
    legend(['COD ',cod_direction],'motion-free frames');
    xlabel('seconds');
    ylabel('mm');
end
disp(['The total number of MFFs are: ',num2str(length(MFF))]);
disp('');
disp('');
disp('% ----------------- %');
disp('% Program Finished! %');
disp('% ----------------- %');


%% Subfunction: "Get_MTPs"
function MTP = Get_MTPs(COD_file_path,start_time,end_time,inflation_factor,n_max,Scout_segments,cod_direction)

    % ----------------------------------------------
    % Additional parameters set up
    % ----------------------------------------------
    n_max_change_scout             = n_max;
    n_scout_segment                = Scout_segments;
    scan_duration_in_use           = end_time-start_time;
    n_max_MFF_use_in_scout_segment = 3;

    allcod                            = dlmread(COD_file_path,' ',1,0);
    if strcmp(cod_direction,'x'); cod = allcod(:,2); end
    if strcmp(cod_direction,'y'); cod = allcod(:,3); end
    if strcmp(cod_direction,'z'); cod = allcod(:,4); end

    cod           = cod(:,1);
    prompt        = allcod(:,end);
    cod_in_use    = cod(start_time+1:end_time,1);
    prompt_in_use = prompt(start_time+1:end_time,1);

    % ----------------------------------------------
    %  Start of SCOUTING: Find the E_tar
    % ----------------------------------------------
    % ----- Plot Figure 1: COD trace for one direction
    figure;
    plot(cod_in_use);
    legend(['COD in ',cod_direction,' direction']); xlabel('seconds'); ylabel('mm');

    % ----- Plot Figure 2: COD trace with initial # of MTPs (n_max_change_scout)
    figure;
    findchangepts(cod_in_use,'MaxNumChanges',n_max_change_scout)
    legend(['COD ',cod_direction],'; model, i.e., mean of COD in each MFF','MTPs (aka changepoints)'); xlabel('seconds'); ylabel('mm');

    [scout_MTP,residue_lowerbound] = findchangepts(cod_in_use,'MaxNumChanges',n_max_change_scout);
    scout_MFF_start                = [1;scout_MTP];
    scout_MFF_end                  = [scout_MTP - 1;scan_duration_in_use];
    scout_MFF_duration             = scout_MFF_end - scout_MFF_start+1;

    disp('scouting motion time points: ');
    for i=1:numel(scout_MTP)
        fprintf(num2str(scout_MTP(i)));fprintf(',');
    end
    disp(['Total residue error from ',num2str(n_max_change_scout),' MTP is: ',num2str(residue_lowerbound)]);

    for idx = 1:numel(scout_MTP)+1
        scout_MFF(idx).start       = scout_MFF_start(idx);
        scout_MFF(idx).end         = scout_MFF_end(idx);
        scout_MFF(idx).duration    = scout_MFF_duration(idx);
        scout_MFF(idx).mean        = mean(cod_in_use(scout_MFF(idx).start:scout_MFF(idx).end));
        scout_MFF(idx).std         = std(cod_in_use(scout_MFF(idx).start:scout_MFF(idx).end));
        scout_MFF(idx).residue     = scout_MFF_duration(idx)*scout_MFF(idx).std.^2;
        scout_MFF(idx).prompt_rate = sum(prompt_in_use(scout_MFF(idx).start:scout_MFF(idx).end))/(scout_MFF(idx).end-scout_MFF(idx).start+1);
    end

    std_scout_MFF_select_in_use = 0;
    E_nm                        = 0;
    duration_scout_segment      = floor(scan_duration_in_use/n_scout_segment);

    % ----------------------------------------------
    % Define scouting segment time mark
    % ----------------------------------------------
    scout_segment_time_mark = [1];
    for idx_scout_segment = 1:n_scout_segment
        scout_segment_time_mark = [scout_segment_time_mark;idx_scout_segment*duration_scout_segment];
    end

    for idx_scout_segment = 1:n_scout_segment % loop over each scouting segment      
        % ----- Find MFF that falls in the current scout_segment
        idx_scout_MFF_select_tmp = intersect(find(scout_MFF_end<scout_segment_time_mark(idx_scout_segment+1)),...
            find(scout_MFF_start>=scout_segment_time_mark(idx_scout_segment)));

        scout_MFF_select_in_use = [];
        if isempty(idx_scout_MFF_select_tmp) % if there is no an complete MFF falls into current scouting segment
            scout_MFF_select_in_use.start       = scout_segment_time_mark(idx_scout_segment);
            scout_MFF_select_in_use.end         = scout_segment_time_mark(idx_scout_segment+1)-1;
            scout_MFF_select_in_use.duration    = duration_scout_segment;
            scout_MFF_select_in_use.mean        = mean(cod_in_use(scout_MFF_select_in_use.start:scout_MFF_select_in_use.end));
            scout_MFF_select_in_use.std         = std(cod_in_use(scout_MFF_select_in_use.start:scout_MFF_select_in_use.end));
            scout_MFF_select_in_use.residue     = scout_MFF_select_in_use.duration*scout_MFF_select_in_use.std.^2;
            scout_MFF_select_in_use.prompt_rate = sum(prompt_in_use(scout_MFF_select_in_use.start:scout_MFF_select_in_use.end))/(scout_MFF_select_in_use.end-scout_MFF_select_in_use.start+1);
        else
            scout_MFF_duration_tmp        = [];
            scout_MFF_select              = scout_MFF(idx_scout_MFF_select_tmp);
            for i                         = 1:numel(scout_MFF_select)
                scout_MFF_duration_tmp(i) = scout_MFF_select(i).duration;
            end

            [~,ascend_idx_scout_MFF_duration] = sort(scout_MFF_duration_tmp);
            n_MFF_use_in_scout_segment_in_use = min([n_max_MFF_use_in_scout_segment,numel(scout_MFF_select)]); % number of MFF in current scout segment will be used
            idx_MFF_scout_in_use              = ascend_idx_scout_MFF_duration(end-n_MFF_use_in_scout_segment_in_use+1:end);
            idx_MFF_scout_in_use              = flip(idx_MFF_scout_in_use); % change to descending

            scout_MFF_select_in_use           = scout_MFF_select(idx_MFF_scout_in_use); % select the MFFs use for scouting
        end

        % ----- Process selected MFF in current scouting segment:
        duration_scout_MFF_select_in_use = 0; % total duration of all MFFs used for scouting in current segment
        residue_scout_MFF_select_in_use  = 0;  % total residue of all MFFs used for scouting in current segment
        for i = 1:numel(scout_MFF_select_in_use)
            residue_scout_MFF_select_in_use  = residue_scout_MFF_select_in_use+scout_MFF_select_in_use(i).residue;
            duration_scout_MFF_select_in_use = duration_scout_MFF_select_in_use+scout_MFF_select_in_use(i).duration;
            std_scout_MFF_select_in_use      = std_scout_MFF_select_in_use + scout_MFF_select_in_use(i).std;
        end
        var_nomotion_predict_cur_segment(idx_scout_segment)=residue_scout_MFF_select_in_use/duration_scout_MFF_select_in_use;
        std_nomotion_predict_cur_segment(idx_scout_segment) = sqrt(var_nomotion_predict_cur_segment(idx_scout_segment));
        residue_nomotion_predict_cur_segment(idx_scout_segment) = duration_scout_segment*var_nomotion_predict_cur_segment(idx_scout_segment); % predicted residue for current segment
        E_nm = E_nm + residue_nomotion_predict_cur_segment(idx_scout_segment);
    end
    E_tar = E_nm * inflation_factor;
    disp(['Predicted total residue for no motion: ',num2str(E_nm)]);
    disp(['Residue target: ',num2str(E_tar)]);fprintf('\n');

    % ----------------------------------------------
    % Find the number of MTP to use (n_tar)
    % ----------------------------------------------
    % ----- Check if E_tar is larger than residue_lowerbound
    if E_tar < residue_lowerbound
        disp('E_tar is found smaller than residue_lowerbound, please raise you MaxNumChanges_ini.');
        disp(['Residue target     : ',num2str(E_tar)]);
        disp(['Residue_lowerbound,: ',num2str(residue_lowerbound)]);fprintf('\n');
        pause;
    else
        for idx_num_change                  = 1:n_max_change_scout
            [~,E_min(idx_num_change)] = findchangepts(cod_in_use,'MaxNumChanges',idx_num_change);
        end
    end

    % ----- Plot Figure 3: Plot MTP vs n
    figure;
    plot(E_min);hold on;
    plot(E_tar*ones(1,n_max_change_scout));hold off;
    title({'minimal total residual error (TRE) for a given number of motion-time-points'});
    legend('residue at different MTPs','residue target');
    xlabel('allowed number of motion time points');
    ylabel('minimal total residual error (TRE)');

    idx_num_change = 1;
    while E_min(idx_num_change) > E_tar
        idx_num_change=idx_num_change+1;
    end
    n_tar = idx_num_change;
    disp(['The number of MTP found after scouting is: ',num2str(n_tar)]);

    % ----- Plot Figure 4: COD trace with new # of MTPs (n_tar, after scouting) 
    figure;
    findchangepts(cod_in_use,'MaxNumChanges',n_tar)
    legend(['COD ',cod_direction],'model, i.e., mean of COD in each MFF','MTPs (aka changepoints)');
    xlabel('seconds');
    ylabel('mm');

    % ----- Finding proper positions of MTPs now that you have n_tar
    [MTP,residue_post_scouting] = findchangepts(cod_in_use,'MaxNumChanges',n_tar);
    disp('The MTPs found after scouting are: ');
    for i = 1:n_tar
        fprintf(num2str(MTP(i)));fprintf(',');
    end
    fprintf('\n');
    disp(['Total residue error, after scouting, from ',num2str(n_tar),' MTP is: ',num2str(residue_post_scouting)]);

    % ----- Calculating motion free frames (MFFs)
    MFF_start    = [1;MTP];
    MFF_end      = [MTP - 1;scan_duration_in_use];
    MFF_duration = MFF_end - MFF_start+1;

    for idx_MFF = 1:numel(MTP)+1
        MFF(idx_MFF).start    = MFF_start(idx_MFF);
        MFF(idx_MFF).end      = MFF_end(idx_MFF);
        MFF(idx_MFF).duration = MFF_duration(idx_MFF);

        % NOTES:
        % - Assign index of scouting segment # to each MFF
        % - If one MFF across two scouting segment, will assign the scouting 
        %   segment number with more MFF duration.
        duration_in_segment=[];
        for idx_scout_segment = 1:n_scout_segment
            duration_in_segment(idx_scout_segment) = numel(intersect((MFF(idx_MFF).start:MFF(idx_MFF).end),...
                scout_segment_time_mark(idx_scout_segment):scout_segment_time_mark(idx_scout_segment+1)));
        end
        [~,belong_idx_scout_segment]   = max(duration_in_segment);
        MFF(idx_MFF).idx_scout_segment = belong_idx_scout_segment;
        MFF(idx_MFF).mean              = mean(cod_in_use(MFF(idx_MFF).start:MFF(idx_MFF).end));
        MFF(idx_MFF).std               = std(cod_in_use(MFF(idx_MFF).start:MFF(idx_MFF).end));
        MFF(idx_MFF).residue           = MFF_duration(idx_MFF)*MFF(idx_MFF).std.^2;
        MFF(idx_MFF).ideal_residue     = MFF_duration(idx_MFF)*var_nomotion_predict_cur_segment(MFF(idx_MFF).idx_scout_segment);
        MFF(idx_MFF).target_residue    = inflation_factor*MFF(idx_MFF).ideal_residue;
        MFF(idx_MFF).prompt_rate       = sum(prompt_in_use(MFF(idx_MFF).start:MFF(idx_MFF).end))/(MFF(idx_MFF).end-MFF(idx_MFF).start+1);
    end
    % ----------------------------------------------
    % End of scouting
    % ----------------------------------------------

    % ----------------------------------------------
    % Calculate summarized results
    % ----------------------------------------------
    % ----- Characterize each MFF
    for idx = 1:numel(MFF)
        score_threshold = [MFF(idx).ideal_residue,(MFF(idx).ideal_residue+MFF(idx).target_residue)/2,MFF(idx).target_residue,2*MFF(idx).target_residue-MFF(idx).ideal_residue];
        if MFF(idx).residue <= score_threshold(1)
            MFF(idx).motionless_score = 1;
        else
            if MFF(idx).residue <= score_threshold(2)
                MFF(idx).motionless_score = 2;
            else
                if MFF(idx).residue <= score_threshold(3)
                    MFF(idx).motionless_score = 3;
                else
                    if MFF(idx).residue <= score_threshold(4)
                        MFF(idx).motionless_score = 4;
                    else
                        MFF(idx).motionless_score = 5;
                    end
                end
            end
        end
    end

    % ----- Summarize
    summary.num_MFF=numel(MFF);
    summary.residue=0;
    summary.ideal_residue = 0;
    summary.target_residue = 0;
    summary.num_MFF_motionless_score1 = 0;
    summary.num_MFF_motionless_score2 = 0;
    summary.num_MFF_motionless_score3 = 0;
    summary.num_MFF_motionless_score4 = 0;
    summary.num_MFF_motionless_score5 = 0;
    summary.weighted_motionless_score = 0;
    tot_weighted_motionless_score = 0;
    tot_duration = 0;

    for idx = 1:numel(MFF)
        summary.residue        = summary.residue + MFF(idx).residue;
        summary.ideal_residue  = summary.ideal_residue + MFF(idx).ideal_residue;
        summary.target_residue = summary.target_residue + MFF(idx).target_residue;

        if MFF(idx).motionless_score==1
            summary.num_MFF_motionless_score1 = summary.num_MFF_motionless_score1+1;
            weighted_motionless_score         = MFF(idx).duration*1;
        end
        if MFF(idx).motionless_score==2
            summary.num_MFF_motionless_score2 = summary.num_MFF_motionless_score2+1;
            weighted_motionless_score         = MFF(idx).duration*2;        
        end
        if MFF(idx).motionless_score==3
            summary.num_MFF_motionless_score3 = summary.num_MFF_motionless_score3+1;
            weighted_motionless_score         = MFF(idx).duration*3;        
        end
        if MFF(idx).motionless_score==4
            summary.num_MFF_motionless_score4 = summary.num_MFF_motionless_score4+1;
            weighted_motionless_score         = MFF(idx).duration*4;        
        end
        if MFF(idx).motionless_score==5
            summary.num_MFF_motionless_score5 = summary.num_MFF_motionless_score5+1;
            weighted_motionless_score         = MFF(idx).duration*5;        
        end
        tot_weighted_motionless_score = tot_weighted_motionless_score+weighted_motionless_score;
        tot_duration                  = tot_duration+MFF(idx).duration;
    end
    summary.weighted_motionless_score = tot_weighted_motionless_score/tot_duration;
    summary.time_preservation         = tot_duration/scan_duration_in_use;

    % ----------------------------------------------
    % Plotting COD trace with MTPs
    % ----------------------------------------------
    MFF_indication = zeros(numel(cod_in_use),1)+min(cod_in_use);
    for idx = 1:numel(MFF)
        MFF_indication(MFF(idx).start:MFF(idx).end-1) = 1*max(cod_in_use);
    end

    % ----- Plot Figure 5: COD trace showing MTPs/MFFs
    figure;
    plot(cod_in_use);hold on;
    plot(MFF_indication);hold on;
    title({['time\_preserve= ',num2str(100*summary.time_preservation),'%  '...
        'Residue\_target= ',num2str(summary.target_residue)],...
        ['MFF\_num = ',int2str(numel(MFF)),...
        ',MaxNumChanges\_ini = ',int2str(n_max_change_scout),...
        ',inflation\_factor = ',num2str(inflation_factor),...
        ',start\_time = ',int2str(start_time)],...
        });
    hold off;
    legend(['COD ',cod_direction],' motion-free frames');
    xlabel('seconds'); ylabel('mm');

    % ----------------------------------------------
    % Saving Residue Information
    % ----------------------------------------------
    Residue.residueEst     = summary.residue;
    Residue.E_tar = summary.target_residue;
    Residue.idealResidue   = summary.ideal_residue;
    save(['COD',cod_direction,'_Residue_Information.mat'],'Residue');
    save(['COD',cod_direction,'_MFF.mat'],'MFF');
    save(['COD',cod_direction,'_MTP.mat'],'MTP');
    close all;

end
    
    