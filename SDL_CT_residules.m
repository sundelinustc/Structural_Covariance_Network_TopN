function SDL_CT_residules(SDL)

% Regress out the confounding factors including age, sex and site
% Cohen's d for between-grou comparisons
% Rank areas showing CT/SA reduction
% Structural covariance matrix based on Pearson's/partial correlation coefficients
% R-to-Z transformation


% if strcmp(SDL.data_type{1}(4:end),'Age') % group x age interaction
%     fprintf('\nFor Group x Age interaction:\n');
%     fprintf('Median_Age=%1.1f yrs\n',median(T.Age));
%     fprintf('N_PTSD_Age>=Median=%d\tN_PTSD_Age<Median=%d\tN_CONT_Age>=Median=%d\tN_CONT_Age<Median=%d\n',...
%         sum(T.Group=='PTSD' & T.Age>=median(T.Age)),...
%         sum(T.Group=='PTSD' & T.Age<median(T.Age)),...
%         sum(T.Group=='CONT' & T.Age>=median(T.Age)),...
%         sum(T.Group=='CONT' & T.Age<median(T.Age)));
% elseif strcmp(SDL.data_type{1}(4:end),'Gender') % group x gender interaction
%     fprintf('\nFor Group x Gender interaction:\n');
%     fprintf('N_PTSD_Male=%d\tN_PTSD_Female=%d\tN_CONT_Male=%d\tN_CONT_Female=%d\n',...
%         sum(T.Group=='PTSD' & T.Gender=='M'),...
%         sum(T.Group=='PTSD' & T.Gender=='F'),...
%         sum(T.Group=='CONT' & T.Gender=='M'),...
%         sum(T.Group=='CONT' & T.Gender=='F'));
% elseif strcmp(SDL.data_type{1}(4:end),'Dep') % group x depression interaction
%     fprintf('\nFor Group x Depression interaction:\n');
%     try
%         fprintf('N_PTSD_Dep=%d\tN_PTSD_NoDep=%d\tN_CONT_Dep=%d\tN_CONT_NoDep=%d\n',...
%             sum(T.Group=='PTSD' & T.Dep==1),...
%             sum(T.Group=='PTSD' & T.Dep==0),...
%             sum(T.Group=='CONT' & T.Dep==1),...
%             sum(T.Group=='CONT' & T.Dep==0));
%     catch
%         fprintf('N_PTSD_Dep=%d\tN_PTSD_NoDep=%d\tN_CONT_Dep=%d\tN_CONT_NoDep=%d\n',...
%             sum(T.Group=='PTSD' & T.Dep=='1'),...
%             sum(T.Group=='PTSD' & T.Dep=='0'),...
%             sum(T.Group=='CONT' & T.Dep=='1'),...
%             sum(T.Group=='CONT' & T.Dep=='0'));
%     end
% else % main effect of group
%     fprintf('\nFor main effect of Group:\n');
%     fprintf('N_PTSD=%d\tN_CONT=%d\n',sum(T.Group=='PTSD'),sum(T.Group=='CONT'));
% end


%% residules through lme regression   
R = []; % matrix containing residuals after lme
k = 0;
for j = 2:149 % per column of data values  
    k = k + 1;
    tbl = T(:,{'Age','Gender','Site'}); tbl.data = T{:,j};
    lme = fitlme(tbl,'data ~ Age + Age^2 + Gender + (1|Site)');
    R(:,k) = residuals(lme);
end
a = T{:,2:149};
figure;
subplot(2,1,1);histogram(a); xlabel('Raw Value'); ylabel('Numbers'); title('Distribution in All Vertices');
subplot(2,1,2);histogram(R); xlabel('Residules After Regressing Out Age, Gender & Site');ylabel('Numbers');
fn = fullfile(SDL.path,'Outputs','CleanData'); 
savefig(fullfile(fn,['Figure Value_Distribution_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
fprintf('Completed: Calculate resudules through regressing out covariates using lme model\n');




% %% Replace cortical thickness outliers (beypnd column mean +/- N*SD) with NaN
% Ncoff = 7; % outliers should beyond mean +/- Ncoff*SD, default is 3, but lead to a lot of outliers in this data
% fprintf('Replacing cortical thickness outliers (beyond column mean +/- %d*SD) with NaN\n',Ncoff);
% for j = col_CT % per column of CT
%     R1 = R(:,j);
%     Rmean = nanmean(R1);
%     Rstd  = nanstd(R1);
%     % look for upper outliers
%     idx1 = []; idx2 = []; idx = [];
%     idx1 = find(R1 > Rmean + Ncoff*Rstd); % index of the upper outliers
%     idx2 = find(R1 < Rmean - Ncoff*Rstd); % index of the lower outliers
%     idx  = [idx1;idx2];
%     
%     if idx % if there is at least an outlier
%         % display the SubjID and cortical area name of the outlier
%         fprintf('--Outliers detected in %s(%d)', T.Properties.VariableNames{j+1},j);
%         for i = 1:length(idx), fprintf('%d(%1.3f<>%1.3f+-%d*%1.3f),',T.SubjID(idx(i)),T{idx(i),j},Tmean,Ncoff,Tstd); end
%         fprintf('\n');
%         % replace outliers
%         T(idx,j) = array2table(NaN);
%     end
% end
% fprintf('Overal values: Min=%1.3f, Mean=%1.3f, Max=%1.3f\n', min(nanmin(T{:,col_CT})), mean(nanmean(T{:,col_CT})), max(nanmax(T{:,col_CT})));






%% Output
fn = fullfile(SDL.path,'Outputs','CleanData'); mkdir(fn);
% if ~strcmp(SDL.ana_type{1},'med') % correlation or regression matrix
%     fn = fullfile(fn,['Data_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% else
    fn = fullfile(fn,['Data_Res_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
% end
save(fn,'T','R','idx1','idx2','d','kk');
% T    --- table containing raw CT values and sbj info
% R    --- residules of CT after lme regression
% idx1 --- row index of PTSD in T & R
% idx2 --- row index of CONT in T & R
% d    --- Cohen's d for PTSD vs CONT
% kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction  
fprintf('Completed: Preprocessed results saved in ->%s\n\n\n',fn);



%%End
end
