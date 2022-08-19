function [data_strct] = makeFTstruct(data_mat,fs,ch_table,sens_info)
% take info from tsv and data matrix and make fieldTrip friendly data
% structure
% 
% data_f...nchan x nsam x n trials data matrix
% 
% fs...sampling frequency
% 
% ch_table...table with fields:
%   name - containing n_chan x 1 cell array with
%            char vectors for channel names
% sens_info...struct containing fields pos and ors, both n_chans x 3 
%             matrices of sensor positions and orientations

data_strct.label = ch_table.name;
data_strct.fsample = fs;
if numel(size(data_mat)) == 2
    Ntrials = 1;
elseif numel(size(data_mat)) == 3
    Ntrials = size(data_mat,3);
else
    error("incorrect number of dimension in data matrix")
end

trial_time = (0:size(data_mat,2)-1)./data_strct.fsample;

for n = 1:Ntrials
    data_strct.trial{n} = data_mat(:,:,n);
    data_strct.time{n} = trial_time;
end

grad.coilpos = sens_info.pos.*100;
grad.coilori = sens_info.ors;
grad.label = data_strct.label;
grad.chanpos = grad.coilpos;
grad.units = 'cm';
data_strct.grad = grad;
end

