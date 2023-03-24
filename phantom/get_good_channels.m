function ch_table = get_good_channels(data_f,ch_table,fs)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Get rid of bad channels
% Get rid of bad channels in data
% ch_table = readtable([path.meg_data,files.channels],'FileType','text','Delimiter','tab');
ch_table.isx = endsWith(ch_table.name,'X');
ch_table.isy = endsWith(ch_table.name,'Y');
ch_table.isz = endsWith(ch_table.name,'Z');
ch_table.slot_no = zeros(height(ch_table),1);
% sanity check
if sum(sum([ch_table.isx,ch_table.isy,ch_table.isz],2)) ~= height(ch_table)
    error('Channel orientation [x,y,z] labels might be wrong!')
end
   
% plot spectra to check for bad channels
[psd_data_all,fxx_d] = get_PSD(data_f,fs);

figure

phs = plot(fxx_d,psd_data_all);set(gca,'YScale','Log');xlim([0,150]);
hold on
for ii = 1:length(phs)
    phs(ii).DisplayName = ch_table.name{ii};
    phs(ii).DataTipTemplate.DataTipRows(end+1)  = dataTipTextRow('Trace',repmat({phs(ii).DisplayName},size(phs(ii).XData))); 
end
plot(fxx_d,mean(psd_data_all,2),'k','LineWidth',2,'DisplayName','Mean')
plotbrowser;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% remove bad channels (data)
% [ch_table] = Bad_Channels(data_f',ch_table,fs);
% %
% %% make sure all bad channels are excluded from noise and MEG data again
% bn = unique(ch_table.name(strcmp(ch_table.status, 'bad')))
% for ii = 1:length(bn) 
%     ch_table.status(find(strcmp(ch_table.name,bn{ii}))) = {'bad'};
%     noise_ch_table.status(find(strcmp(noise_ch_table.name,bn{ii}))) = {'bad'};
% end
% 
% %% Check final matrices are the same
% 
% % writetable(ch_table,[path.meg_data,filename,files.channels,'_new'],...
% %         'WriteRowNames',true,'Delimiter','tab','FileType','text')
% % 
% % ch_table = readtable([path.meg_data,filename,files.channels,'_new'],...
% %         'Delimiter','tab','FileType','text');
% %     
% % writetable(noise_ch_table,sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),...
% %         'WriteRowNames',true,'Delimiter','tab','FileType','text')
% % noise_ch_table = readtable(sprintf('%s%s_%s%s_new',path.noise,files.noise,ses,files.noise_channels),...
% %         'Delimiter','tab','FileType','text');
% % disp("Removing bad channels")
% bad_chans_data = [find(startsWith(ch_table.status,'bad'))];
% % bad_chans_noise = [find(startsWith(noise_ch_table.status,'bad'))];
% 
% ch_table_ = ch_table;
% data_f_ = data_f;
% % noise_ch_table_ = noise_ch_table;
% % noise_data_f_ = noise_data_f;
% 
% 
% ch_table_(bad_chans_data,:) = [];
% data_f_(bad_chans_data,:) = [];
% % noise_ch_table_(bad_chans_noise,:) = [];%clear noise_ch_table
% % noise_data_f_(bad_chans_noise,:) = [];
% % 
% % if size(data_f_,1) ~= size(noise_data_f_,1)
% %     error('Channel counts not the same!')
% % end
% % 
% % if ~all([ch_table_.name{:}] == [noise_ch_table_.name{:}])
% %     error('Noise and data channels not in same order!')
% % end
% %%
% [psd_data_all,fxx_d] = get_PSD(data_f_,fs);
% 
% % [psd_data_all_n,fxx_n] = get_PSD(noise_data_f_,fs);
% 
% figure
% % phs = plot(fxx_n,psd_data_all_n,'r');set(gca,'YScale','Log');xlim([0,150]);hold on
% % for ii = 1:length(phs);phs(ii).DisplayName = noise_ch_table.name{ii};end
% % hold on
% phs = plot(fxx_d,psd_data_all);set(gca,'YScale','Log');xlim([0,150]);
% for ii = 1:length(phs);phs(ii).DisplayName = ch_table.name{ii};end
% plotbrowser;