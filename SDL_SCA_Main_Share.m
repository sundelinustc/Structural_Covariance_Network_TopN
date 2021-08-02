
%==============================================
% This script is for top-n structural covariance analyses (SCA)
% by Dr. Delin Sun, in Morey Lab, Duke University, 07-26-2019
%==============================================
clear;
close all;

%% Parameters

cd ..
SDL.path = pwd; % path of the project
cd('Scripts');
SDL.raw = fullfile(SDL.path,'Original','data_example.xlsx'); % raw data file
SDL.out = fullfile(SDL.path,'Outputs_ComBat_Reduction');mkdir(SDL.out); % find the top-n regions of PTSD < Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_Increase');mkdir(SDL.out); % find the top-n regions of PTSD > Non-PTSD
% SDL.out = fullfile(SDL.path,'Outputs_ComBat_NoDiff');mkdir(SDL.out); % find the top-n regions of PTSD ~ Non-PTSD

% data type & analyses type
% CT = cortical thickness, SA = surface area
% corr = correlation, partialcorr = partial correlation,
% regr = regression, regrcov = regression with covariates (time consuming)
% med = mediation, medcov = mediation with covariates (time consuming)
% und = undirected connection, dir = directed connection
% ---data_type, ana_type, connection with direction or not, predefined X, Y and M
Ana = {
    'CT', 'corr',        'und', {},{},{},'';
    'SA', 'corr',        'und', {},{},{},'';
    };

SDL.group = {'PTSD','CONT'};
SDL.N = 5000; % number of permutation, 5000 for publication purpose
SDL.Nr =300; % number of searching matched solutions, 300 for publicaiton purpose

disp(SDL);


%% Analyses
for i = 1:size(Ana,1) % per pre-defined analysis
    SDL.data_type = Ana(i,1); % CT or SA
    SDL.ana_type  = Ana(i,2); % corr, partialcorr or med
    SDL.XYM       = Ana(i,4:7); % for mediation analyses only
    
    
    SDL_Load_Clean(SDL); % load and clean raw data
    SDL_ComBat(SDL); % Harmonization data across study sites using ComBat
    SDL_LME_residules(SDL); % residules after LME to exclude age, age^2, sex, and mean CT or SA
    SDL_CTSA_Comp(SDL); % find the top-n regions of PTSD < Non-PTSD
    
    SDL_SCA_TopN(SDL); % SCA withn the top-N areas, Notice: time consuming for 5k permutation!!
    
end