function SDL_ComBat(SDL)

% ComBat to remove the effects of study sites


%% add path of combat scripts
addpath(fullfile(SDL.path,'Scripts','ComBat'));

%% load original data
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Data_Clean_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
load(fn); fprintf('Loaded: %s\n',fn);

% if strcmp(SDL.data_type{1},'CT_Dep') || strcmp(SDL.data_type{1},'SA_Dep')
%     T = T(~isnan(T.Dep),:); % must have depression data
% elseif strcmp(SDL.data_type{1},'CT_PTSDsev10') || strcmp(SDL.data_type{1},'SA_PTSDsev10')
%     T = T(~isnan(T.CAPS) & T.CAPStype==4,:); % must have CAPS-IV data
% else
% end

%% ComBat protocal
batch = T.Site; % site
dat = T{:,2:149}'; % values, Num_voxel x Num-sbj
age = T.Age;
T.Sex(strcmp(T.Gender,'M')) = 1; T.Sex(strcmp(T.Gender,'F')) = 2;
sex = dummyvar(T.Sex); %sex, 1=male, 2=female, dummyvar only accepts positive integer
dep = dummyvar(T.Dep+1);
caps = T.CAPS;
T.Grp(strcmp(T.Group,'PTSD')) = 2; T.Grp(strcmp(T.Group,'CONT')) = 1;
disease = dummyvar(T.Grp); % diagnosis, 1=CONT, 2=PTSD, dummyvar only accepts positive integer

[FILEPATH,NAME,EXT] = fileparts(SDL.out);
if strcmp(NAME,'Outputs_ComBat')
%     if strcmp(SDL.data_type{1},'CT_Dep') || strcmp(SDL.data_type{1},'SA_Dep')
%         mod = [age sex(:,2) dep(:,2) disease(:,2)];
%     elseif strcmp(SDL.data_type{1},'CT_PTSDsev10') || strcmp(SDL.data_type{1},'SA_PTSDsev10')
%         mod = [age sex(:,2) caps];
%     else % for main, age inter, sex inter
        mod = [age sex(:,2) disease(:,2)];
%     end
elseif strcmp(NAME,'Outputs')
    mod = disease(:,2);
else
    mod = [];
end
% disease 1st column is control, 2nd is patients
% sex 1st column is male, 2nd is female
% dep 2nd column is subjects with depression

% data_harmonized = combat(dat, batch, mod, 0); % non-parametric adjustments
data_harmonized = combat(dat, batch, mod, 1); % parametric adjustments

% plot figure
get(0,'Factory'); set(0,'defaultfigurecolor',[1 1 1]); % change background and marginal areas color into white
figure;
subplot(1,2,1);histogram(T{:,2:149}); xlabel(SDL.data_type{1}); ylabel('Numbers'); title('Cleaned Data');
subplot(1,2,2);histogram(data_harmonized); xlabel(SDL.data_type{1});ylabel('Numbers');title('Harmonized Data Using ComBat');
fdir = fullfile(SDL.out,SDL.data_type{1}); 
savefig(fullfile(fdir,['Figure Value_Cleaned_vs_ComBat_',SDL.data_type{1},'_',SDL.ana_type{1},'.fig']));
fprintf('Completed: ComBat to harmonize data across study sites\n');

% replace the raw values in the table
T{:,2:149} = data_harmonized';

%% save data
fdir = fullfile(SDL.out,SDL.data_type{1});
fn = fullfile(fdir,['Data_ComBat_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
save(fn,'T');
fprintf('\nSaved: ComBat data --> %s\n', fn);

%% remove path of combat scripts
rmpath(fullfile(SDL.path,'Scripts','ComBat'));


%% End
end