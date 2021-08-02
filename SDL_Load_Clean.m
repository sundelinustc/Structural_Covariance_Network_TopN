function SDL_Load_Clean(SDL)
% Load and clean data
% load raw data (from Rakesh), and clean data by replacing 0, NaN and outliers

%% Load raw data
fname = SDL.raw;
T = readtable(fname,'sheet',SDL.data_type{1}(1:2));
if strcmp(SDL.out(end-4:end),'splt1')
    T = T(1:2:size(T,1),:);
elseif strcmp(SDL.out(end-4:end),'splt2')
    T = T(2:2:size(T,1),:);
end
T = T(T.Included==1,:); % remove the red rows
T = T(T.Gender==1 | T.Gender==2,:); % keep only clear males & females
fprintf('Completed: Load<- %s\n', SDL.raw);

%% Rename Groups
T.Grp(T.Group==1) = {'PTSD'};
T.Grp(T.Group==2) = {'CONT'};
T.Group = T.Grp; T.Grp = [];
fprintf('Completed: Rename PTSD & CONT\n');

%% Re-organize data table into ['CONT';'PTSD']
T1 = T(strcmp('PTSD',T.Group),:); 
T2 = T(strcmp('CONT',T.Group),:);
T = [T1;T2];
fprintf('Completed: Re-organize data table into PTSD and CONT\n');

%% Rename Gender
T.Gen(T.Gender==1) = {'M'};
T.Gen(T.Gender==2) = {'F'};
T.Gender = T.Gen; T.Gen = [];
T = T(strcmp('M',T.Gender) | strcmp('F',T.Gender),:); % exclude the subjects missing gender values
fprintf('Completed: Re-organize values of Gender\n');

%% Replace cortical thickness 0 with NaN
col_CT     = 2:149; % the columns containing cortical thickness data 
T1         = T{:,col_CT};
T1(T1==0)  = NaN; % the cortical thickness data are in columns 2-149
T(:,col_CT)= array2table(T1);
fprintf('Completed: Replace all cortical values of 0s with NaN\n');


%% Replace NaN with Column Mean per group
T1 = T(strcmp('PTSD',T.Group),:);
T2 = T(strcmp('CONT',T.Group),:);
for j = [col_CT,150] % per column of CT/SA (2~149), Age (150)
    T1(isnan(T1{:,j}),j) = array2table(nanmean(T1{:,j})); 
    T2(isnan(T2{:,j}),j) = array2table(nanmean(T2{:,j})); 
end
T = [T1;T2];
fprintf('Completed: Replace cortical thickness of NaN with Column Mean per group\n');

%% Some Covariates set as Categorical variable (for lme)
% T.Group  = nominal(T.Group);
% T.Gender = nominal(T.Gender);
% T.Site   = nominal(T.Site);
T.Dep    = T.Dep_yes_no_1_Yes_0_no; 
T.PTSDsev = T.PTSD_HighVersusLow_1_high_0_low;
% try 
%     sum(T.Dep == 1);
% catch
%     T.Dep = nominal(T.Dep);
% end


%% Output
fdir = fullfile(SDL.out,SDL.data_type{1}); mkdir(fdir);
% if ~strcmp(SDL.ana_type{1},'med') % correlation or regression matrix
%     fn = fullfile(fdir,['Data_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% else
    fn = fullfile(fdir,['Data_Clean_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
% end
save(fn,'T');
% T    --- table containing raw CT values and sbj info
fprintf('Saved: Cleaned data saved in ->%s\n\n\n',fn);

%% End
end