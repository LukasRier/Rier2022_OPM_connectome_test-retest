clear all
close all
clc

% set project dir to directory containing data from https://doi.org/10.5072/zenodo.1134455

project_dir =  'path/to/project/'

%% Generate virtual electrodes and source power for MEG and empty room noise data
for sub_i =1:10
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        resting_VE_noise_proj_func(sub,ses,project_dir)
    end
end
%% Estimate functional connectivity
for sub_i = 1:10
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        resting_conn_func(sub,ses,project_dir)
    end
end
