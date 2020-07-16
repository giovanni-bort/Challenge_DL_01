 %function save_CINC_ECG_IMG
 
% salva un file IMG  per paziente in CINC20_ECG_IMG1
% nome file: ECG_Xnnnn_ddd.mat   X:AQISHE    nnnn: 0001    ddd:0167  diagnosi(es. 01+67)
%
 
ECG_mean=[ I_mean II_mean III_mean V1_mean V2_mean V3_mean V4_mean V5_mean V6_mean ];
ECG_TMP=[ I;II;III; V1; V2; V3; V4; V5; V6];
if(size(ECG_TMP,2)<5000), ECG_TMP(1,5000)=0;end
ECG_PRO=ECG_TMP(1:9,1:5000)';
%ECG_CINC.dati=ECG_PRO(:);
% STRUCT(num_file)
ECG_CINC.dati=[ECG_PRO(:)' ECG_mean];
 ECG_CINC.diagn=STRUCT(num_file).diagn;
 ECG_CINC.ind_diagn=22;
ECG_CINC.Fs=Hz;
[ind_dia_star,TEMP2,TEMP3]=check_DIAGN_star(STRUCT(num_file).diagn);

% % [KNEW_NAME,K_TYPE,K_NUM]=extract_info_from_name(H_recording);
% % if(ne(K_NUM,num_file)),fprintf('***** DIVERSI***');end
% % fprintf('NEW:%s  K_NUM=%6.0f  num_file:%6.0f  type: %s  %s \n',KNEW_NAME,K_NUM,num_file,K_TYPE,file_key);
% fprintf('size DIA_*:%6.0f%6.0f\n',size(ind_dia_star));
%  file_00X=fullfile(ECG_DL_directory,['ECG_',file_key,num2str(num_file,'%05.0f_'),num2str(STRUCT(num_file).ind_diagn,'%02.0f')]);
%  file_00X=fullfile(['ECG_',file_key,num2str(num_file,'%05.0f_'),num2str(STRUCT(num_file).ind_diagn,'%02.0f')]);

%file_00X=fullfile(['ECG_',file_key,num2str(num_file,'%05.0f_'),num2str(ind_dia_star,'%02.0f')]); % ERROR**
file_00X=fullfile(['ECG_',file_key,num2str(num_file,'%05.0f_'),num2str(ind_dia_star,'%02.0f') ]);

 
 %----------------- OPT_IMG=1  --------------------------------------------
 opt_IMG=1;
if(opt_IMG==1)
   imgLoc =ECG_DL_directory;
% %    if ~exist(imgLoc, 'dir')
% %         mkdir(imgLoc)
% %    end
   fprintf('Save ECG images in : %s (Wavelet scalogram)\n',imgLoc);
    data=ECG_CINC.dati(:);
    Fs=ECG_CINC.Fs;
    [~,signalLength] = size(data);
    signalLength = numel(data);
    fprintf(' -> %6.0f samples ',signalLength);
    fb = cwtfilterbank('SignalLength',signalLength, 'SamplingFrequency',Fs,'VoicesPerOctave',12);

    cfs = abs(fb.wt(data));
    im = ind2rgb(im2uint8(rescale(cfs)),jet(128));
    
    %imgLoc = fullfile(imageRoot,char(labels(ii)));
%     imgLoc = imageRoot;
    
%     imFileName = strcat(num2str(ii,'%04.0f_'),char(ALL_LABELS(ii)),'_',num2str(ii),'.jpg');        % ARR_23.jpg
%     imFileName = strcat(num2str(ii,'%04.0f_'),num2str(ii),'.jpg');  
%     imFileName = strrep(files(ii).name,'.mat','.jpg');
    opt_QUAL=75;
    imFileName = [file_00X '.jpg'];
    imwrite(imresize(im,[224 224]),fullfile(imgLoc,imFileName),'quality',opt_QUAL);
    fprintf('IMG:   %s ',imFileName); fprintf(' %s ',ECG_CINC.diagn{:}); fprintf(' %3.0f ',ECG_CINC.ind_diagn); fprintf('\n');

end
   