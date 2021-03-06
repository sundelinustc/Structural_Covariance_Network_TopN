function SDL_Euler(SDL)

% calculation of Euler number
fname = fullfile(SDL.path,'Original','euler_stats.xlsx');
T = readtable(fname,'sheet','Euler');
T = T(:,1:5); % name (site), orig_id, Group, Euler-lh, Euler-rh
Sites = unique(T.Name); % 28 sites here

% left hemi
fprintf('\nLeft Hemisphere\n');
for i = 1:size(Sites,1) % per site
    SiteName = Sites{i};
    T0 = T(strcmp(T.Name,SiteName),:); % data of this site
    T1 = T0(T0.Group==1,'Euler_lh'); T1(any(ismissing(T1),2), :) = []; % PTSD, remove NaN
    T2 = T0(T0.Group==2,'Euler_lh'); T2(any(ismissing(T2),2), :) = []; % CONT, remove NaN
    if ~size(T1,1) | ~size(T2,1) % if there is NO data
        fprintf('Site,PTSD_mean,PTSD_std,CONT_mean,CONT_std,t,df,p:\t%s\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n',SiteName);
    else
        [H,p,CI,stats] = ttest2(T1.Euler_lh,T2.Euler_lh);
        fprintf('Site,PTSD_mean,PTSD_std,CONT_mean,CONT_std,t,df,p:\t%s\t%d\t%d\t%d\t%d\t%1.3f\t%d\t%1.3f\n',...
            SiteName,mean(T1.Euler_lh),std(T1.Euler_lh),mean(T2.Euler_lh),std(T2.Euler_lh),...
            stats.tstat,stats.df,p);
    end
end


% right hemi
fprintf('\nRight Hemisphere\n');
for i = 1:size(Sites,1) % per site
    SiteName = Sites{i};
    T0 = T(strcmp(T.Name,SiteName),:); % data of this site
    T1 = T0(T0.Group==1,'Euler_rh'); T1(any(ismissing(T1),2), :) = []; % PTSD, remove NaN
    T2 = T0(T0.Group==2,'Euler_rh'); T2(any(ismissing(T2),2), :) = []; % CONT, remove NaN
    if ~size(T1,1) | ~size(T2,1) % if there is NO data
        fprintf('Site,PTSD_mean,PTSD_std,CONT_mean,CONT_std,t,df,p:\t%s\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\tNaN\n',SiteName);
    else
        [H,p,CI,stats] = ttest2(T1.Euler_rh,T2.Euler_rh);
        fprintf('Site,PTSD_mean,PTSD_std,CONT_mean,CONT_std,t,df,p:\t%s\t%d\t%d\t%d\t%d\t%1.3f\t%d\t%1.3f\n',...
            SiteName,mean(T1.Euler_rh),std(T1.Euler_rh),mean(T2.Euler_rh),std(T2.Euler_rh),...
            stats.tstat,stats.df,p);
    end
end

%% End
end