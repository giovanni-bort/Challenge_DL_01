function [list_ok,files_ok,full_files_ok,List_Diagn]=select_train_records(dir_IMG,n_MAX_rec,n_DIA_chk,NORM_chk,opt_prn)
%-------INPUT:---------
% 1 dir_IMG   : directory of IMG
% 2 n_MAX_rec : numero max di record per classe diagnostico
% 3 n_DIA_chk :
% 4 NOrm_CHK  : 
% 5 opt_prn   : 0:no_print    1:short    2:long 
%----OUTPUT:--------------
% 1 list_ok       : lista numero di files selezionati di dir_IMG
% 2 files_ok      : nomi files selezionati
% 3 full_files_ok : full namees dei files selezionati
% 4 List_Diagn    : diagnosi dei files selezionati
%

if(nargin<2),  n_MAX_rec=0; end
if(nargin<3),  n_DIA_chk=0; end
if(nargin<4),  NORM_chk=0;  end       %NORM_chk : 1: tolgo NORMAL se presenti altre diagnosi
if(nargin<5),  opt_prn=0;   end

DIAGN_star=[270492004 	164889003 	164890007 	426627000 	713427006 	713426002 	445118002   39732003 	164909002 ...
            251146004 	698252002 	10370003 	284470004 	427172004  164947007 	111975006 	164917005 	47665007  ...
            59118001 	427393009 	426177001   426783006 	427084000 	63593006 	164934002 	59931005 	17338001 9999999 ];
HDIAGN_star ={ 'IAVB' ;   'AF'    ; 'AFL'     ; 'Barady'  ; 'CRBBB'    ; 'IRBBB'  ; 'LAnFB'   ; 'LAD'     ; 'LBBB' ; ...
               'LQRSV' ; 'NSIVCB' ; 'PR'      ;  'PAC'    ;  'PVC'     ;  'LQT'   ;  'LPR'    ; 'QaB'     ;  'RAD' ; ...
               'RBBB'  ;  'SA'    ; 'SB'      ; 'SNR'     ; 'STach'    ; 'SVPB'   ; 'TAb'     ; 'Tinv'    ; 'VEB' ; '*Other*' } ;
N28=numel(DIAGN_star);


IMG_files=dir([dir_IMG filesep 'ECG_*.jpg']);
n_files=numel(IMG_files);
if(opt_prn>0), fprintf(' presenti %6.0f files\n',n_files);end

    II=[1:10 1050:1060  7000:7010  10000:10010 15000:15010  18000:18010  20000:20010   22000:22010  24000:24010 25000:25010];
Tot_Dia=0;Tot_Dia(N28)=0;Cross_Dia=0;Cross_Dia(N28,N28)=0;DB_Dia=0;DB_Dia(6,N28)=0;
Ref_Dia(1).n_files=[];Ref_Dia(N28).n_files=[];
%-----Loop n_files --------------------------------
tot_summ=0;List_Diagn=[];
for i_f=1:n_files
    [N_rec,N_DB,typ_DB,diagn]=extract_ndiagn(IMG_files(i_f).name);
    
    if(NORM_chk>0)
        if(numel(diagn)>1), diagn=setdiff(diagn,22);end     % tolgo NORMAL in presenza di più diagnosi
    end
    
    tot_summ=tot_summ+numel(diagn);
    
    
    for j_d=diagn
       Ref_Dia(j_d).n_files(end+1)=i_f;
    end
    Tot_Dia(diagn)=Tot_Dia(diagn)+1;
    
    DB_Dia(N_DB,diagn)=DB_Dia(N_DB,diagn)+1;
    
    Cross_Dia(diagn,diagn)=Cross_Dia(diagn,diagn)+1;

    III=intersect(i_f,II);
    List_Diagn{i_f}=diagn;
    if(numel(III)>0 & opt_prn>1)
        fprintf('%5.0f ',i_f);
      fprintf('%-25s ',IMG_files(i_f).name);
      fprintf('N:%6.0f DB:%1.0f Key:%s',N_rec,N_DB,typ_DB);
      fprintf(' diagn: '); fprintf('%3.0f',diagn);
      fprintf('\n');
    end
end
%-------- end Loop n_files ----------------------------------------
fprintf(' presenti: %6.0f files  , %6.0f diagnosi\n',n_files,tot_summ);

if(opt_prn>0)
 if(NORM_chk==0),fprintf('---------------------------------------------------------\n');end
 if(NORM_chk>0) ,fprintf('----------------- *** NORMAL check ***-------------------------\n');end
 fprintf('    N_tot   SNOMED_code      code    ');fprintf('  DB_%1.0f',1:6);fprintf('\n')
 for i_cl=1:N28
    fprintf('%3.0f  %4.0f',i_cl,Tot_Dia(i_cl));
    fprintf('%13.0f',DIAGN_star(i_cl));
    fprintf('   %3.0f ',i_cl);
    fprintf('%-8s',HDIAGN_star{i_cl});
    fprintf('%6.0f',DB_Dia(:,i_cl));
    fprintf('\n');
 end
 fprintf('tot%6.0f                            ',sum(Tot_Dia)); fprintf('%6.0f',sum(DB_Dia,2)); 
 fprintf('  -> %6.0f',sum(sum(DB_Dia)));    fprintf('\n');
 fprintf('------without class 28:"other"---------------------------------------------------\n');
  fprintf('tot%6.0f                            ',sum(Tot_Dia(1:N28-1))); fprintf('%6.0f',sum(DB_Dia(:,1:N28-1),2)); 
 fprintf('  -> %6.0f',sum(sum(DB_Dia(:,1:N28-1))));    fprintf('\n');
 fprintf('---------------------------------------------------------\n');
fprintf('  ');fprintf('%5.0f',1:N28);fprintf('\n');   
 for i_cl=1:N28
    fprintf('%2.0f',i_cl);
    fprintf('%5.0f',Cross_Dia(i_cl,:));
    fprintf('\n');
 end
 fprintf('---------------------------------------------------------\n');
end

if(n_MAX_rec==0),
   n_MAX_rec=input('numero massimo records x classe [100] ');
   if(isempty(n_MAX_rec)),n_MAX_rec=100;end
end
list_tot=[];
classi_chk=1:N28-1;
if(n_DIA_chk>0),classi_chk=n_DIA_chk;end
for i_cl=classi_chk
    NN=numel(Ref_Dia(i_cl).n_files);
    if(NN>0)
       rng('default');
       new_ind=randperm(NN,min(n_MAX_rec,NN));
       Ref_Dia(i_cl).new_files=Ref_Dia(i_cl).n_files(new_ind);
       list_tot=[list_tot Ref_Dia(i_cl).new_files];
    end
    if(NN==1),fprintf('%6.0f   NN=%6.0f\n',i_cl,NN);end 
end


list_ok=unique(list_tot);

ii_0=find(list_tot==0);if(numel(ii_0)>0),fprintf('ERROR *** Zero:');fprintf('%6.0f',ii_0);fprintf('\n'); end

files_ok=[];full_files_ok=[];
for i_f=list_ok
    files_ok{end+1}=IMG_files(i_f).name;
     full_files_ok{end+1}=[IMG_files(i_f).folder '\' IMG_files(i_f).name];

end

List_Diagn=List_Diagn(list_ok);
N_tot_DIA=0;N_tot_DIA(28)=0;     K_tot_DIA=0;K_tot_DIA(28)=0;
for i=1:numel(list_ok)
    N_tot_DIA(List_Diagn{i})=N_tot_DIA(List_Diagn{i})+1;
     K_tot_DIA(List_Diagn{i})=K_tot_DIA(List_Diagn{i})+1;
     if(numel(List_Diagn{i})>1 & numel(intersect(List_Diagn{i},22))>0),K_tot_DIA(22)=K_tot_DIA(22)-1;end
end


if(opt_prn>0)
 for i_cl=classi_chk
    fprintf('%3.0f',i_cl);
    fprintf('%6.0f',numel(Ref_Dia(i_cl).n_files));
    fprintf(' ->  %6.0f',numel(Ref_Dia(i_cl).new_files));
    fprintf('  -->%6.0f',N_tot_DIA(i_cl));
    fprintf('  corr: %6.0f',K_tot_DIA(i_cl));
    
    nn_zero=numel(find(Ref_Dia(i_cl).new_files==0));
    if(nn_zero>0),fprintf(' ZERO:%6.0f',nn_zero);end
    fprintf('\n');
 end
end

fprintf('** (N_max=%4.0f) Selected total files=%6.0f    ',n_MAX_rec,numel(list_tot));
fprintf('  uniques:%6.0f  n_dia=%6.0f    corr:%6.0f\n',numel(list_ok),sum(N_tot_DIA),sum(K_tot_DIA));

end
%----------------------------------------------------------------------
%===================================================================================================
%----------------------------------------------------------------------
function [N_rec,N_DB,typ_DB,diagn]=extract_ndiagn(stringa)
% ECG_XNNNN_DDDD.jpg
ID1=0;ID2=0;ID3=0;ID4=0;
N_rec=str2num(stringa(6:10));
typ_DB=stringa(5);
N_DB=0;
if(typ_DB=='A'),N_DB=1;end
if(typ_DB=='Q'),N_DB=2;end
if(typ_DB=='I'),N_DB=3;end
if(typ_DB=='S'),N_DB=4;end
if(typ_DB=='H'),N_DB=5;end
if(typ_DB=='E'),N_DB=6;end
if(N_DB==0),fprintf('***** ERROR ****\n');end
k=0;
for i=12:2:length(stringa)-4
    k=k+1;
    diagn(k)=str2num(stringa(i:i+1));
end


end
