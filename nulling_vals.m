clearvars

project_dir = '/path/to/data/folder';
if ~exist(project_dir,'dir')
    error('Set project directory!')
end
fc_dir = ['data',filesep,'derivatives',filesep,'field_control',filesep];
results_dir = [project_dir,'results',filesep];

subs = ['001';'001';'002';'002';'003';'003';'004';'004';'005';'005';'006';'006';...
    '007';'007';'009';'009';'010';'010'];
% Load each file with model coefficients
for fi = 1:length(subs)
   run = 2-mod(fi,2);
   L = load(sprintf('%s%ssub-%s\\sub-%s_model_coefficients_null_%d.mat',project_dir,fc_dir,subs(fi,:),subs(fi,:),run));
   
   % Root sum square magnitude of the mapped field
   null_B_mags(fi) = L.B_mag;
   % Root sum square magnitude of the mapped field gradients
   null_grad_mags(fi) = sqrt(sum(L.coeffs(5:8).^2));
end

% separate values into runs
run1_Bmags = null_B_mags(1:2:18);
run2_Bmags = null_B_mags(2:2:18);
run1_grad_mags = null_grad_mags(1:2:18);
run2_grad_mags = null_grad_mags(2:2:18);

% Calculate means and standard deviations
mean_Bmags_run1 = mean(run1_Bmags)
mean_Bmags_run2 = mean(run2_Bmags)

std_Bmags_run1 = std(run1_Bmags)
std_Bmags_run2 = std(run2_Bmags)

mean_grad_mags_run1 = mean(run1_grad_mags)
mean_grad_mags_run2 = mean(run2_grad_mags)

std_grad_mags_run1 = std(run1_grad_mags)
std_grad_mags_run2 = std(run2_grad_mags)

% Save to file in results directory
f = fopen([results_dir,'Nulling_results.txt'],'w');
fprintf(f,'%s\n           Field mapping results\n%s\n',repmat('-',1,45),repmat('-',1,45));
fprintf(f,'Null 1:\nRSS field magnitude = %1.4f +/- %1.4f nT\n',mean_Bmags_run1,std_Bmags_run1);
fprintf(f,'RSS field gradients = %1.4f +/- %1.4f nT/m\n',mean_grad_mags_run1,std_grad_mags_run1);
fprintf(f,'Null 2:\nRSS field magnitude = %1.4f +/- %1.4f nT\n',mean_Bmags_run2,std_Bmags_run2);
fprintf(f,'RSS field gradients = %1.4f +/- %1.4f nT/m\n',mean_grad_mags_run2,std_grad_mags_run2);
fclose(f);

fprintf('%s\nField mapping results\n%s\n',repmat('-',1,45),repmat('-',1,45));
fprintf('Null 1:\nRSS field magnitude = %1.4f +/- %1.4f nT\n',mean_Bmags_run1,std_Bmags_run1);
fprintf('RSS field gradients = %1.4f +/- %1.4f nT/m\n',mean_grad_mags_run1,std_grad_mags_run1);
fprintf('Null 2:\nRSS field magnitude = %1.4f +/- %1.4f nT\n',mean_Bmags_run2,std_Bmags_run2);
fprintf('RSS field gradients = %1.4f +/- %1.4f nT/m\n',mean_grad_mags_run2,std_grad_mags_run2);
