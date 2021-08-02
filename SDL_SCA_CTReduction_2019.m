
%==============================================
% This script is for structural covariance analyses (SCA) on cortical thickness (CT) & surface area (SA) changes,
% by Dr. Delin Sun, in Morey Lab, Duke University, 07-26-2019
%==============================================
clear;
close all;

%% Parameters
% SDL.path = '/Volumes/dusom_morey/Data/Lab/Delin/Projects/ENIGMA_SCA_CT_Reduction'; %  my Mac desktopathp in the lab
% SDL.path = '\\duhsnas-pri.dhe.duke.edu\dusom_morey\Data\Lab\Delin\Projects\ENIGMA_SCA_CT_Reduction'; % path my Lenova laptop at home
% SDL.BoxPath = '/Users/ds366/Documents/MATLAB/SDLToolBox'; % path for most of the toolboxes in the lap desktop
% SDL.BoxPath = 'C:\Users\sunde\Documents\MATLAB\SDLToolBox'; % path for most of the toolboxes at home

cd ..
SDL.path = pwd; % path of the project
cd('Scripts');
SDL.raw = fullfile(SDL.path,'Original','PGC_Data_Covariates_check_for_delin_withCAPS.xlsx'); % raw data file

SDL.out = fullfile(SDL.path,'Outputs_ComBat_Reduction');mkdir(SDL.out); % find the top-n regions of PTSD < Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_Increase');mkdir(SDL.out); % find the top-n regions of PTSD > Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_NoDiff');mkdir(SDL.out); % find the top-n regions of PTSD ~ Non-PTSD

% SDL.out = fullfile(SDL.path,'Outputs_ComBat_Reduction_splt2');mkdir(SDL.out); % find the top-n regions of PTSD < Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_Increase');mkdir(SDL.out); % find the top-n regions of PTSD > Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_NoDiff');mkdir(SDL.out); % find the top-n regions of PTSD ~ Non-PTSD

SDL.Boxpath = fullfile('C:\Users\sunde\Documents\MATLAB','SDLToolBox');

% data type & analyses type
% CT = cortical thickness, SA = surface area
% corr = correlation, partialcorr = partial correlation, 
% regr = regression, regrcov = regression with covariates (time consuming)
% med = mediation, medcov = mediation with covariates (time consuming)
% und = undirected connection, dir = directed connection
% ---data_type, ana_type, connection with direction or not, predefined X, Y and M
Ana = {
%     'CT', 'corr',        'und', {},{},{},''; 
%     'SA', 'corr',        'und', {},{},{},'';

%     'CT_PTSDsev10', 'corr',    'und', {},{},{},'';
%     'SA_PTSDsev10', 'corr',    'und', {},{},{},'';

%     'CT_Dep', 'corr',    'und', {},{},{},'';
%     'SA_Dep', 'corr',    'und', {},{},{},'';

%     'CT_Gender', 'corr', 'und', {},{},{},''; 
%     'SA_Gender', 'corr', 'und', {},{},{},'';
 
%     'CT_Age10', 'corr',    'und', {},{},{},''; % arbitarily divide into 6 groups 0~20, 20~30 ...
%     'SA_Age10', 'corr',    'und', {},{},{},'';
    };

SDL.group = {'PTSD','CONT'};
SDL.N = 5000; % number of permutation, 5000 for publication purpose
SDL.Nr =300; % number of searching matched solutions, 300 for publicaiton purpose

disp(SDL);


%% Analyses
% for i = 1:size(Ana,1) % per pre-defined analysis
%     SDL.data_type = Ana(i,1); % CT or SA
%     SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
%     SDL.XYM       = Ana(i,4:7); % for mediation analyses only
%     
%     
%         SDL_Load_Clean(SDL); % load and clean raw data
%         SDL_ComBat(SDL); % Harmonization data across study sites using ComBat
%         SDL_LME_residules(SDL); % residules after LME to exclude age, age^2, sex, and mean CT or SA
%         SDL_CTSA_Comp(SDL); % find the top-n regions of PTSD <, > or = Non-PTSD   
%     
%     if strcmp(SDL.data_type{1},'CT') || strcmp(SDL.data_type{1},'SA')
%         SDL_SCA_All(SDL); % SCA within the whole brain, Notice: time consuming for 5k permutation!!
%         SDL_SCA_TopN(SDL); % SCA withn the top-N areas, Notice: time consuming for 5k permutation!!
%     elseif strcmp(SDL.data_type{1},'CT_PTSDsev10') || strcmp(SDL.data_type{1},'SA_PTSDsev10')
%         % %     SDL_SCA_All_PTSDsev10(SDL);
%         SDL_SCA_TopN_PTSDsev10(SDL);
%     elseif strcmp(SDL.data_type{1},'CT_Age10') || strcmp(SDL.data_type{1},'SA_Age10')
%         %         SDL_SCA_All10(SDL);  % SCA within the whole brain for multiple groups basd on age
%         SDL_SCA_TopN10(SDL); % SCA withn the top-N areas for multiple groups basd on age
%         %     SDL_power(SDL); % power analysis for Michael Debellis's grant
%     else
%         %         SDL_SCA_All4(SDL);  % SCA within the whole brain for 4 groups
%         SDL_SCA_TopN4(SDL); % SCA withn the top-N areas for 4 groups
%     end
%     
% end

% SDL_Final_Fig1(SDL);
% SDL_Final_Fig2(SDL);
% SDL_Final_Fig3(SDL);
% SDL_Final_Fig4(SDL);
% SDL_Final_Fig5(SDL);
% SDL_Final_Fig6(SDL);
% SDL_Final_Fig7(SDL);


% Table2(SDL); % PTSD vs. random, non-PTSD vs. random, and PTSD vs. non-PTSD
% Table3(SDL); % PTSD severity (CAPS scores)
% Table4(SDL); % PTSD x depression
% Table5(SDL); % PTSD x gender
% Table6(SDL); % PTSD modified by Age

% SDL_Euler(SDL); % calculation of Euler number
SDL_SCA_LeaveNSiteOut(SDL); % leave N site(s) out analyses