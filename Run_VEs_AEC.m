% clearvars
close all
clc
% Conn_analysis_MJB('001','001');
% Conn_analysis_MJB('001','002');
% clear all
% Conn_analysis_MJB('002','001');
% Conn_analysis_MJB('002','002');
% clear all
% Conn_analysis_MJB('003','001');
% Conn_analysis_MJB('003','002');
% clear all
% Conn_analysis_MJB('004','001');
% Conn_analysis_MJB('004','002');
% clear all
% Conn_analysis_MJB('005','001');
% Conn_analysis_MJB('005','002');
% clear all
% Conn_analysis_MJB('006','001');
% Conn_analysis_MJB('006','002');
% clear all
% Conn_analysis_MJB('007','001');
% Conn_analysis_MJB('007','002');
% clear all
% Conn_analysis_MJB('008','001');
% Conn_analysis_MJB('008','002');
% clear all
% Conn_analysis_MJB('009','001');
% Conn_analysis_MJB('009','002');
% clear all
% Conn_analysis_MJB('010','001');
% Conn_analysis_MJB('010','002');


% restoredefaultpath
% ft_defaults

for sub_i = 1:10
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        data_cleaning_test(sub,ses)
    end
end
