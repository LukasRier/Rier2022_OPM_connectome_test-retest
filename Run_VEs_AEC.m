% clearvars
close all
clc

% restoredefaultpath
% ft_defaults

for sub_i = 1:10
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        resting_conn_func(sub,ses)
    end
end
