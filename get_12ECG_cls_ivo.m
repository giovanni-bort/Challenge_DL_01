
function [scores,out_labels]=get_12ECG_cls_ivo(ECG,Header,model,classes)
        
%addpath(genpath('My_Tools/'))
        

out_labels=0;out_labels(1:numel(classes))=0;
scores=out_labels;

[H_recording,Total_time,H_num_leads,H_Fs,H_gain,H_age,H_sex]=extract_data_from_header(Header);

    fprintf('File,time');fprintf(' %s %8.2f ',H_recording,Total_time);
    fprintf('leads:%6.0f Fs:%8.1f age:sex%6.0f%6.0f\n',H_num_leads,H_Fs,H_age,H_sex);
    fprintf('Gain:');fprintf('%8.1f',H_gain);
    fprintf(' max:%8.0f',max(max(ECG)));
    fprintf('\n');

global out_labels_1 out_labels_2 KK_ERROR
     try
%    Hz=500;
   Hz=H_Fs;
   Freq=Hz; 
    fprintf('Version LAST_ECG_07_09_function -- size(ECG)=%6.0f%8.0f\n',size(ECG));

    j=1;
     filename=' ';
   
%     do_04_21_function
%     LAST_ECG_07_09_function    
    last_05_27_function    

 valuta_trained_NN
    
    DIAGN_star=[270492004 	164889003 	164890007 	426627000 	713427006 	713426002 	445118002   39732003 	164909002 ...
            251146004 	698252002 	10370003 	284470004 	427172004  164947007 	111975006 	164917005 	47665007  ...
            59118001 	427393009 	426177001   426783006 	427084000 	63593006 	164934002 	59931005 	17338001 9999999 ];
classes_dbl=0;
        for i_cl=1:numel(classes)
    classes_dbl(i_cl)=str2double(classes{i_cl});
end
[IC1,IC2,IC3]=intersect(classes_dbl,DIAGN_star);

my_label_ext=0; my_label_ext(28)=0;
my_scores_ext=0; my_scores_ext(28)=0;
my_label_ext([1:4 6:23 25:26])=my_label;
my_scores_ext([1:4 6:23 25:26])=my_scores;

%  for i=[1:4 6:23 25:26],k=k+1;classes_star{k}=num2str(i,'%02.0f');end
fprintf('label:');fprintf('%6.1f',my_label);fprintf('\n');
fprintf('labEX:');fprintf('%6.1f',my_label_ext);fprintf('\n');
fprintf('score:');fprintf('%6.3f',my_scores);fprintf('\n');
fprintf('scoEX:');fprintf('%6.3f',my_scores_ext);fprintf('\n');

out_labels=0;out_labels(1:numel(classes))=0;
scores=0;    scores(1:numel(classes))=0;
if(numel(IC1)>0)
  out_labels(IC2)=my_label_ext(IC3);
  scores  = my_scores_ext(IC3);
else
    fprintf('***** non common labels ***** classes: %6.0f\n',numel(classes));
end
    
% % % out_labels = [AF(j) , I_AVB(j) , LBBB(j) , Normal(j) , PAC(j) , PVC(j), RBBB(j),  STD(j), STE(j), ];
% % % 
% % % scores=out_labels/(max(sum(out_labels),1) ) ;
% % % 

     catch MSG_ERROR
        num_file=1;file_key='A';
paz_error=sprintf('ERROR: file:%6.0f Key:%s; ',num_file,file_key);
text_error=[paz_error ' ' MSG_ERROR.message ' ' MSG_ERROR.stack(1).name ' line:' num2str(MSG_ERROR.stack(1).line)];
   KK_ERROR=KK_ERROR+1;
   fprintf('------ ERORR n. %6.0f -------------\n',KK_ERROR);
   fprintf(':ERORR: %s\n',text_error);

save_annotation(text_error,'err_annotations.txt');
    end



end

% 
% function plot(a,b,c,d,e,f,g,h,i,j,k,l,m,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)
% end
% function subplot(a,b,c,d)
% end
% function figure(a,b,c)
% end
% function title(a,b,c,d)
% end
% function clf
% end
% function hold (a,b,c)
% end
% function axis(a,b,c,d,e,f)
% end
% function view(a,b,c,d,e,f,g,h)
% end
% function xlabel(a,b,c,d,f,g)
% end
% function ylabel(a,b,c,d,e,f,g)
% end
