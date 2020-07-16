
function [tmp_data]=train_12ECG_cls_ivo(ECG,Header,STRUCT,ECG_DL_directory,num_file,file_key)

%addpath(genpath('My_Tools/'))

tmp_data=0;
scores=0;
out_labels=0;
[H_recording,Total_time,H_num_leads,H_Fs,H_gain,H_age,H_sex]=extract_data_from_header(Header);

    fprintf('File,time');fprintf(' %s %8.2f ',H_recording,Total_time);
    fprintf('leads:%6.0f Fs:%8.1f age:sex%6.0f%6.0f\n',H_num_leads,H_Fs,H_age,H_sex);
    fprintf('Gain:');fprintf('%8.1f',H_gain);
    fprintf(' max:%8.0f',max(max(ECG)));
    fprintf('\n');

     %Version do_04_02  copiede 11.4.20
global out_labels_1 out_labels_2 KK_ERROR
     try
%    Hz=500;
   Hz=H_Fs;
    
    fprintf('Version last_05_27_function -- size(ECG)=%6.0f%8.0f\n',size(ECG));

    j=1;
     filename=' ';
   
%     do_04_21_function
    last_05_27_function    

%fprintf('file name:%s\n',name);

 % {'AF'} {'I-AVB'} {'LBBB'} {'Normal'} {'PAC'} {'PVC'} {'RBBB'} {'STD'} {'STE'}
 
% % % out_labels = [AF(j) , I_AVB(j) , LBBB(j) , Normal(j) , PAC(j) , PVC(j), RBBB(j),  STD(j), STE(j), ];
% % % 
% % % scores=out_labels/(max(sum(out_labels),1) ) ;
% % % 
% % % fprintf('do-labels:');fprintf('%8.0f',out_labels); fprintf('\n');
% % % fprintf('do-scores:');fprintf('%8.3f',scores); fprintf('\n');
save_CINC_ECG_IMG
% out_labels_1=out_labels_1 + out_labels;

% ******  For Deep learning  ***** DO NOT REMOVE ****
%     K_OPT_NN=1;
%    valuta_trained_NN;     
%     fprintf('NN1-labels:');fprintf('%8.0f',out_labels); fprintf('\n');
%     fprintf('NN1-scores:');fprintf('%8.3f',scores); fprintf('\n');
  
% % %         K_OPT_NN=2;
% % %    valuta_trained_NN;     
% % %     fprintf('NN2-labels:');fprintf('%8.0f',out_labels); fprintf('\n');
% % %     fprintf('NN2-scores:');fprintf('%8.3f',scores); fprintf('\n');
% %     K_OPT_NN=3;
% %    valuta_trained_NN;     
% %     fprintf('NN3-labels:');fprintf('%8.0f',out_labels); fprintf('\n');
% %     fprintf('NN3-scores:');fprintf('%8.3f',scores); fprintf('\n');


% out_labels_2=out_labels_2 + out_labels;
%     fprintf('tot_labels1-1:');fprintf('%8.0f',out_labels_1); fprintf('\n');
%     fprintf('tot_labels1-2:');fprintf('%8.0f',out_labels_2); fprintf('\n');

     catch MSG_ERROR
        
paz_error=sprintf('ERROR: file:%6.0f Key:%s; ',num_file,file_key);
text_error=[paz_error ' ' MSG_ERROR.message ' ' MSG_ERROR.stack(1).name ' line:' num2str(MSG_ERROR.stack(1).line)];
   KK_ERROR=KK_ERROR+1;
   fprintf('------ ERORR n. %6.0f -------------\n',KK_ERROR);
%    fprintf(':ERORR: %s\n',text_error);

save_annotation(text_error,'err_annotations.txt');
    end



end

