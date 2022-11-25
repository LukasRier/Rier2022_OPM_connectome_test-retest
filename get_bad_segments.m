function [artifacts] = get_bad_segments(data_mat,thresh_val)
Ntrials = size(data_mat,3);
std_data_mat = squeeze(std(data_mat,[],2));
std_data_mat_cent = std_data_mat - repmat(mean(std_data_mat,2),1,size(std_data_mat,2));
std_trials = std(std_data_mat_cent,[],2);
threestd_trials_mat = repmat(thresh_val.*std_trials,1,size(std_data_mat,2));
logic_mat = std_data_mat_cent > threestd_trials_mat;
bad_trials_vec = sum(logic_mat,1);
bad_trials = find(bad_trials_vec > 1);
artifacts = [];
trl_samps = size(data_mat,2);
for ind = 1:length(bad_trials)
    artifacts(ind,1) = 1000+(bad_trials(ind)-1).*trl_samps;
end
artifacts = [artifacts,artifacts+trl_samps-2000];
end