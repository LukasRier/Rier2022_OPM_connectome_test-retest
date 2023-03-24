function [Data] = read_cMEG_data(filename)
            % Basic read in of raw data with Channel names and session
            % info. Does NOT apply conversion factor. This should be
            % recorded by the operator. Gain x1 has conversion 2.7 V/nT,
            % adjust accordingly for other gain settings (e.g. x0.33 Gain =
            % 0.9 V/nT)
            if strcmpi(num2str(filename(end-4:end)),'.cMEG')
                filename = filename(1:end-5);
            end

            disp('Get data')
            fname = [filename '.cMEG'];
            fid = fopen(fname,'rb','ieee-be');
            finfo = dir(fname);
            fsize = finfo.bytes;
            Adim_conv = [2^32; 2^16; 2^8; 1]; % Array dimension conversion table

            disp('Preallocation')
            I = 0;
            while fsize ~= ftell(fid)
%                     disp(I)
                dims = [];
                for n = 1:2 % 2 dimensions (Nch x Time)
                    temp = fread(fid,4);
                    temp = sum(Adim_conv.*temp);
                    dims = [dims,temp];
                end
                I = I + 1;
                temp = fread(fid,prod(dims),'double',0,'ieee-be');  % Skip the actual data (much faster for some reason)
            end
            fseek(fid,0,'bof'); % Reset the cursor
            data1 = repmat({NaN*ones(dims)},I,1);  % Preallocate space, assumes each section is the same

            disp('Load and parse data')
            for j = 1:I
                dims = [];  % Actual array dimensions
                for n = 1:2 % 2 dimensions (Nch x Time)
                    temp = fread(fid,4);
                    temp = sum(Adim_conv.*temp);
                    dims = [dims,temp];
                end
                clear temp n
                temp = fread(fid,prod(dims),'double',0,'ieee-be');

                % Reshape the data into the correct array configuration
                for n = 1:2
                    temp = reshape(temp,dims(2-n+1),[])';
                end
                data1{j} = temp;  % Save the data
            end
            fclose(fid);  % Close the file

            disp('Reshape into sensible order')
            data = NaN*zeros(size(data1{1},2),size(data1{1},1).*size(data1,1));
            for n = 1:size(data1,1)
                clc
                %     disp(100*n/length(data1))
                data_range = [(n-1)*size(data1{1},1)+1:n*size(data1{1},1)];
                data(:,data_range) = [data1{n}']; % data in Nch+triggers+1 x time, but channel one is the time
            end

            disp('Read session info')
            Session_info_file = [filename(1:end-4) '_SessionInfo.txt'];
            fidSI = fopen(Session_info_file);
            finfoSI = dir(Session_info_file);
            fsizeSI = finfoSI.bytes;
            count = 0;
            while fsizeSI ~= ftell(fidSI)
                count = count + 1;
                Session_info{count,1} = fgetl(fidSI);
            end
            fclose(fidSI);

            disp('Get channel names')
            for n = 1:size(Session_info,1)
                if strfind(Session_info{n},'Sensor Names:')
                    Chan_info_string_id = n;
                end
            end
            Chan_info_string = Session_info{Chan_info_string_id};
            cn1 = strsplit(Chan_info_string,',');
            cn11 = strsplit(cn1{1},': ');
            Chan_names = [cn11(2) cn1(2:end)]';

            disp('Create QZFM_data structure')
            Data.Session_info = Session_info;
            Data.samp_frequency = round(1/(data(1,2)-data(1,1)));
            Data.nsamples = size(data,2);
            Data.time = linspace(0,size(data,2)./Data.samp_frequency,size(data,2));
            Data.elapsed_time = Data.nsamples./Data.samp_frequency;
            Data.Chan_names = Chan_names;
            Data.data = data(2:end,:);
        end