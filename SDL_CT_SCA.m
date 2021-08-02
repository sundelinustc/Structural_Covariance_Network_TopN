
function SDL_CT_SCA(SDL)


% Mean SCA values of top-N areas with CT/SA reduction
% Mean SCA values of 5k random sets of N areas

%% Loading information
fn = fullfile(SDL.path,'Outputs','CleanData'); 
% if ~strcmp(SDL.ana_type{1},'med') % correlation or regression matrix
%     fn = fullfile(fn,['Data_',SDL.data_type{1},'_',SDL.ana_type{1},'.mat']);
% else
    fn = fullfile(fn,['Data_Res_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
% end
load(fn);
% T    --- table containing raw CT values and sbj info
% R    --- residules of CT after lme regression
% idx1 --- row index of PTSD in T & R
% idx2 --- row index of CONT in T & R
% d    --- Cohen's d for PTSD vs CONT
% kk   --- the list of ranked area No. for PTSD vs CONT based on d, i.e. kk(1) showes the area with the largest CT reduction  
fprintf('Completed: Load preprocessed results<- %s\n\n\n',fn);


%% Structural Covariance Matrix
if strcmp(SDL.ana_type{1},'corr')
    func = @corr;
elseif strcmp(SDL.ana_type{1},'partialcorr')
    func = @partialcorr;
else
end
S0 = SDL_SCA(R,idx1,idx2,func,0);     % SC of original PTSD & CONT data
S1 = SDL_SCA(R,idx1,idx2,func,SDL.N);  % SC of N times' permutations of PTSD & CONT labels, for whole-brain betwen-group contrast
fn = fullfile(SDL.path,'Outputs','CleanData',['Data_SC_',SDL.data_type{1},'_',SDL.ana_type{1},'_',SDL.XYM{4},'.mat']);
S0.M1 = [];S0.M2 = []; S1.M1 = [];S1.M2 = []; % to save space cause Matlab has difficulty to save data larger than 2GB
save(fn,'S0','S1','kk');
fprintf('Completed: Structural covariance matrix saved in -> %s\n\n\n',fn);


%% End
end

