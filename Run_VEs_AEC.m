clear all
close all
clc


% project_dir = pwd;project_dir(end-6:end)=[];
project_dir =  '/net/cador/data_local/Lukas/movie/'
% for sub_i =10
%     sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
%     for ses_i = 1:2
%         ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
%         resting_VE_noise_proj_func(sub,ses,project_dir)
%     end
% end
%%
for sub_i = 5:6
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        resting_conn_func(sub,ses,project_dir)
    end
end

% for sub_i =1:10
%     sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
%     for ses_i = 1:2
%         ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
%         resting_VE_pow(sub,ses,project_dir)
%     end
% end