function [classes, STRUCT]=driver_train(input_directory, ECG_DL_directory,KK_ini,KK_end)

    if(nargin<3),KK_ini=1; end
    if(nargin<4),KK_end=0; end
	% Find files.
    input_files = {};
    for f = dir(input_directory)'
        if exist(fullfile(input_directory, f.name), 'file') == 2 && f.name(1) ~= '.' && all(f.name(end - 2 : end) == 'mat')
            input_files{end + 1} = f.name;
        end
    end

    
    % read number of unique classes
    [classes,STRUCT] = get_classes(input_directory,input_files);
    for ii=1:numel(STRUCT)
        [~,ind_diagn]=ismember(STRUCT(ii).diagn,classes);
        STRUCT(ii).ind_diagn=ind_diagn;
    end

    if(KK_ini>0)
       if ~exist(ECG_DL_directory, 'dir')
        mkdir(ECG_DL_directory)
    end

 

    % Load model.
     disp('Loading 12ECG model...')
     model = load_12ECG_model();
%model=STRUCT;
    % Iterate over files.
    disp('Predicting 12ECG labels...')
    num_files = length(input_files);
    if(KK_end==0),KK_end=num_files;end
    for i = KK_ini:KK_end  %% num_files
 %   for i = 1:num_files
         file_key=input_files{i}(1:1);

        disp(['    ', num2str(i), '/', num2str(num_files), '...' , input_files{i},  ' Key: ', file_key])

        % Load data.
        file_tmp=strsplit(input_files{i},'.');
        tmp_input_file = fullfile(input_directory, file_tmp{1});
        [data,header_data] = load_challenge_data(tmp_input_file);
        [current_score,current_label] = run_12ECG_classifier(data,header_data,classes,model,STRUCT,ECG_DL_directory,i,file_key);

%         save_challenge_predictions(output_directory,file_tmp{1}, current_score, current_label,classes);
	
    end


    disp('Done.')
    end   % ---- if KK_ini>0 ---
    
end


function [data,tlines] = load_challenge_data(filename)

        % Opening header file
        fid=fopen([filename '.hea']);
        if (fid<=0)
                disp(['error in opening file ' filename]);
        end

        tline = fgetl(fid);
        tlines = cell(0,1);
        while ischar(tline)
            tlines{end+1,1} = tline;
            tline = fgetl(fid);
        end
        fclose(fid);

        f=load([filename '.mat']);
        try
                data = f.val;
        catch ex
                rethrow(ex);
        end

end

% find unique number of classes
function [classes,STRUCT] = get_classes(input_directory,files)
	
	classes={};
    STRUCT=[];
	num_files = length(files);
	k=1;
    	for i = 1:num_files
            STRUCT(i).num=i;
            STRUCT(i).file=files{i};
            
            
		g = strrep(files{i},'.mat','.hea');
		input_file = fullfile(input_directory, g);
	        fid=fopen(input_file);
	        tline = fgetl(fid);
        	tlines = cell(0,1);

		while ischar(tline)
        	    tlines{end+1,1} = tline;
	            tline = fgetl(fid);
			if startsWith(tline,'#Dx')
				tmp = strsplit(tline,': ');
				tmp_c = strsplit(tmp{2},',');
				for j=1:length(tmp_c)
                                        STRUCT(i).diagn{j}=tmp_c{j};

		                	idx2 = find(strcmp(classes,tmp_c{j}));
		                	if isempty(idx2)
                	        		classes{k}=tmp_c{j};
                        			k=k+1;
                			end
				end
			break
        		end
		end
        	fclose(fid);
	end
	classes=sort(classes)
end



function save_challenge_predictions(output_directory,recording, scores, labels,classes)

	output_file = ([output_directory filesep recording '.csv']);

	Total_classes = strjoin(classes,','); %insert commaas
	%write header to file
	fid = fopen(output_file,'w');
	fprintf(fid,'#%s\n',recording);
	fprintf(fid,'%s\n',Total_classes);
	fclose(fid);

	%write data to end of file
	dlmwrite(output_file,labels,'delimiter',',','-append','precision',4);
        dlmwrite(output_file,scores,'delimiter',',','-append','precision',4);

end

