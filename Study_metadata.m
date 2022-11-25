clear all
close all
clc
all_sub = {};
all_ses = {};
N_data_ch = [];
N_bad_ch = [];
raw_data_dur = [];
Noise_dur = [];
N_bad_epoch = [];
N_bad_ica_comps = [];
project_dir = pwd;project_dir(end-6:end)=[];
project_dir = '/net/cador/data_local/Lukas/movie/';
for sub_i =1:10
    sub = sprintf('%3d',sub_i);sub(sub == ' ') = '0'
    for ses_i = 1:2
        ses = sprintf('%3d',ses_i);ses(ses == ' ') = '0'
        all_sub = cat(1,all_sub,sub);
        all_ses = cat(1,all_ses,ses);
        addpath([project_dir,'scripts',filesep,'fieldtrip-20190212'])
        addpath([project_dir,'scripts'])
        addpath([project_dir,'scripts',filesep,'Beamformer',filesep,''])
        ft_defaults;
        datadir = [project_dir,'data',filesep];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        exp_type = 'task-movie';
        filename = ['sub-',sub,'_',exp_type,'_','run-',ses];
        path.main = [datadir,'sub-',sub,filesep];
        path.ICA = [datadir,'derivatives',filesep,'ICA',filesep,'sub-',sub,filesep];
        path.cleaning = [datadir,'derivatives',filesep,'cleaning',filesep,'sub-',sub,filesep];
        path.meshes = [datadir,'derivatives',filesep,'sourcespace',filesep,'sub-',sub,filesep];
        path.VEs = [datadir,'derivatives',filesep,'VEs',filesep,'sub-',sub,filesep];
        path.noiseVEs = [datadir,'derivatives',filesep,'noiseVEs',filesep,'sub-',sub,filesep];
        path.pow = [datadir,'derivatives',filesep,'POW',filesep,'sub-',sub,filesep];
        
        if ~exist(path.ICA,'dir'); mkdir(path.ICA);end
        if ~exist(path.cleaning,'dir'); mkdir(path.cleaning);end
        if ~exist(path.VEs,'dir'); mkdir(path.VEs);end
        if ~exist(path.noiseVEs,'dir'); mkdir(path.noiseVEs);end
        if ~exist(path.pow,'dir'); mkdir(path.pow);end
        
        path.data = [path.main,'meg',filesep];
        path.noise = [path.main,'noise',filesep];
        path.mri = [path.main,'anat',filesep];
        files.mri = ['sub-',sub,'.nii'];
        files.meshes = ['sub-',sub,'_meshes.mat'];
        files.VEs = ['sub-',sub,'_','run-',ses,'_VE'];
        files.noise = ['sub-',sub,'_task-noise'];
        files.noiseVEs = ['sub-',sub,'_','run-',ses,'_noiseVE'];
        files.pow = ['sub-',sub,'_','run-',ses,'_POW'];
        files.voxlox = ['sub-',sub,'_voxlox.mat'];
        files.channels = '_channels.tsv';
        files.noise_channels = ['_channels.tsv'];
        files.sens_order = '_sensor_order.mat';
        files.ICA = [path.ICA,filename];
        files.cleaning = [path.cleaning,filename];
        S.mri_file = [path.meshes,'sub-',sub,'_',filesep,'x.nii']; % only need path for beamformer and single shell function. bit hacky, FIX
        % global files
        cd(path.data)
        load([path.data,filename,'_meg.mat'],'fs','data','triggers')
        load([datadir,'derivatives',filesep,'helmet',filesep,'M1p5_adult.mat'])
        noisedata = load(sprintf('%s%s.mat',path.noise,files.noise));
        noise_fs = noisedata.fs; noise_data = noisedata.data; clear noisedata
        % basefilename = [path.main,'ses-',ses,'/Matt/'];
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %% Preproc
        % Mean correct
        data = data - mean(data,1);
        noise_data = noise_data - mean(noise_data,1);
        N_data_ch = cat(1,N_data_ch,size(data,2));

        % Notch filter
        for harms = [50,100,150,200,120]
            Wo = harms/(fs/2);  BW = Wo/35;
            [b,a] = iirnotch(Wo,BW);
            if fs ~= noise_fs; error("Noise and data not sampled at the same frequency!");end
            disp('Applying Notch filter')
            data = filter(b,a,data,[],1);
            noise_data = filter(b,a,noise_data,[],1);
        end
        
        % bandpass filter for viewing
        disp('Applying 1-150 Hz bandpass filter')
        
        hp = 1;
        lp = 150;
        [b,a] = butter(4,2*[hp lp]/fs);
        data_f = [filtfilt(b,a,data)]';
        noise_data_f = [filtfilt(b,a,noise_data)]';
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Get rid of bad channels
        %% Get rid of bad channels in data
        
        if ~exist([path.data,filename,files.channels,'_new'],'file')
            error("Use 'get_good_channels.m for this dataset first")
        else
            ch_table = readtable([path.data,filename,files.channels,'_new'],...
                'Delimiter','tab','FileType','text');
        end
        
        % Bad channels in noise recording
        if ~exist(sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),'file')
            error("Use 'get_good_channels.m for this dataset first")
            
        else
            noise_ch_table = readtable(sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),...
                'Delimiter','tab','FileType','text');
        end
        %%
        
        load([path.data,filename,files.sens_order],'T')
        for sl_i = 1:height(T)
            if ~isempty(T{sl_i,1}{1})
                ch_table.slot_no(startsWith(ch_table.name,T{sl_i,1}{1}(1:2))) = sl_i;
            end
        end
        %%
        
        %%
        % remove from data matrix
        disp("Removing bad channels")
        bad_chans_data = [find(startsWith(ch_table.status,'bad'))];
        bad_chans_noise = [find(startsWith(noise_ch_table.status,'bad'))];
        N_bad_ch = cat(1,N_bad_ch,size(bad_chans_data,1));

        ch_table_ = ch_table;
        data_f_ = data_f;
        noise_ch_table_ = noise_ch_table;
        noise_data_f_ = noise_data_f;
        
        
        ch_table_(bad_chans_data,:) = [];
        data_f_(bad_chans_data,:) = [];
        noise_ch_table_(bad_chans_noise,:) = [];%clear noise_ch_table
        noise_data_f_(bad_chans_noise,:) = [];
        
        if size(data_f_,1) ~= size(noise_data_f_,1)
            error('Channel counts not the same!')
        elseif ~all([ch_table_.name{:}] == [noise_ch_table_.name{:}])
            error('Noise and data channels not in same order!')
        else
            ch_table = ch_table_;
            data_f = data_f_;
            noise_ch_table = noise_ch_table_;
            noise_data_f = noise_data_f_;
            clear ch_table_ data_f_ noise_ch_table_ noise_data_f_
        end
                        
        %% Epoch data
        % segment using trigger
        disp("Chopping data into 5s epochs")
        tr = diff(triggers(:,1)>2);
        start_sample1 = find(tr == 1);
        start_sample = start_sample1(1);
        end_sample = find(tr == -1);
        duration = floor((end_sample - start_sample)/fs);
        
        % data = data(start_sample+1:start_sample+duration*fs,:);
        data_f = data_f(:,start_sample+1:start_sample+duration*fs);
        raw_data_dur = cat(1,raw_data_dur,duration);
        Noise_dur = cat(1,Noise_dur,size(noise_data,1)/fs);
        
        %% Check for bad epochs
        
        % chop into short segments
        epoch_length = 5;
        data_f_mat = reshape(data_f,[size(data_f,1),epoch_length*fs,duration/epoch_length]);
        clear data_f
        load([files.cleaning,'_vis_artfcts.mat'],'vis_artfcts')
        
        % automatic artifact rejection
        thresh_val = 3;
        auto_artfcts = get_bad_segments(data_f_mat,thresh_val);
        fprintf('Found %d artifacts using a threshold of %d std. deviations.\n',...
            size(auto_artfcts,1),thresh_val);
        
        all_artfct = [vis_artfcts;auto_artfcts];
        
        N_bad_epoch = cat(1,N_bad_epoch,size(all_artfct,1));
        
        clear data_f_mat
        
        %% ICA
        
        disp("Loading bad coponents, topographies and old unmixing matrix")
        load([files.ICA,'_bad_ICA_comps.mat'],'bad_comps')
        N_bad_ica_comps = cat(1,N_bad_ica_comps,length(bad_comps));
    end
end
study_meta_data = table(all_sub,all_ses,N_data_ch,N_bad_ch,raw_data_dur,...
    Noise_dur,N_bad_epoch,N_bad_ica_comps,N_data_ch-N_bad_ch,...
    100.*(N_data_ch-N_bad_ch)./N_data_ch,raw_data_dur-(N_bad_epoch.*epoch_length),'VariableNames',...
    {'sub','ses','N_data_ch','N_bad_ch','raw_data_dur',...
    'Noise_dur','N_bad_epoch','N_bad_ica_comps','n_good_chans','good_chans_prct','clean_data'});
writetable(study_meta_data,[project_dir,'Study_meta_data.xlsx'])