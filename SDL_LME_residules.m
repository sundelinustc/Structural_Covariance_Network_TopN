function SDL_LME_residules(SDL)

% residules through lme regression   



%% load ComBat data
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Data_ComBat_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn);
fprintf('\nLoaded: ComBat data <-- %s\n', fn);

%% Linear Mixed Model (LME)
R = []; % matrix containing residuals after lme
T1 = T;
T.Gender = nominal(T.Gender);
T.mean = mean(T{:,2:149},2); % mean across areas

% % between-group difference of mean values
% V_PTSD = mean(T{strcmp(T.Group,'PTSD'),2:149},2);
% V_CONT = mean(T{strcmp(T.Group,'CONT'),2:149},2);
% [H,P,CI,STATS] = ttest2(V_PTSD,V_CONT);

for j = 1:148 % per column of data values
    tbl = T(:,{'Age','Gender','mean'});
    tbl.data = T{:,j+1}; % data are in column 2-149
    if strcmp(SDL.data_type{1},'CT_Age10') || strcmp(SDL.data_type{1},'SA_Age10')
        lme = fitlme(tbl,'data ~ Gender + mean');
    elseif strcmp(SDL.data_type{1},'CT_Gender') || strcmp(SDL.data_type{1},'SA_Gender')
        lme = fitlme(tbl,'data ~ Age + Age^2 + mean');
    else
        lme = fitlme(tbl,'data ~ Age + Age^2 + Gender + mean');
    end
    R(:,j) = residuals(lme);
end
T = T1;

%% plot
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure('DefaultAxesFontSize',18)
subplot(1,2,1);histogram(T{:,2:149}); xlabel(SDL.data_type{1}); ylabel('Numbers'); title('ComBat Data');
subplot(1,2,2);histogram(R); xlabel(SDL.data_type{1});ylabel('Numbers');title(['Residules After Regressing Out Age, Age^2, Gender & mean ',SDL.data_type{1}]);
fdir = fullfile(SDL.out,SDL.data_type{1}); 
savefig(fullfile(fdir,['Figure Residules_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
fprintf('Completed: Calculate resudules through regressing out covariates using lme model\n');

% replace the raw values in the table
T{:,2:149} = R; % residules

%% Save data
fdir = fullfile(SDL.out,SDL.data_type{1}); 
fn = fullfile(fdir,['Data_Residuals_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
save(fn,'T');
fprintf('\nSaved: Residuals --> %s\n', fn);


%% End
end