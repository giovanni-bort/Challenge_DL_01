function [IA,IB,IC]= check_DIAGN_star(diagn_to_chk)
%
% IA : nuovo numero diagnosi tra 1:28    (28=other)
% IB : abbreviazioni diagnosi star
% IC : codici SNOMED diagnosi star
%

%    **************************** NON MODIFICARE ORDINE DIAGN_star ====> codici diagnosi per training CNN ************************     
%    **************************  (original in check_diagn_star) ------------------------------------------------------------------
DIAGN_star=[270492004 	164889003 	164890007 	426627000 	713427006 	713426002 	445118002   39732003 	164909002 ...
            251146004 	698252002 	10370003 	284470004 	427172004  164947007 	111975006 	164917005 	47665007  ...
            59118001 	427393009 	426177001   426783006 	427084000 	63593006 	164934002 	59931005 	17338001 9999999 ];
HDIAGN_star ={ 'IAVB' ;   'AF'    ; 'AFL'     ; 'Barady'  ; 'CRBBB'    ; 'IRBBB'  ; 'LAnFB'   ; 'LAD'     ; 'LBBB' ; ...
               'LQRSV' ; 'NSIVCB' ; 'PR'      ;  'PAC'    ;  'PVC'     ;  'LPR'    ;  'LQT'   ; 'QaB'     ;  'RAD' ; ...
               'RBBB'  ;  'SA'    ; 'SB'      ; 'SNR'     ; 'STach'    ; 'SVPB'   ; 'TAb'     ; 'Tinv'    ; 'VEB' ; '*Other*' } ;


%removed:             x                                                        x       x 
%                     A                       B  C              A              B       C 
%       [ 1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 ];
IA_fin= [ 1  2  3  4 19  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 13 25 26 14 28 ];

 A=diagn_to_chk;B=0;for i=1:numel(A), B(i)=str2double(A{i});end
[C,IA]=intersect(DIAGN_star,B);
if(numel(C)==0),C=999999; IA=28;end
fprintf(': ');fprintf('%3.0f',IA); 
dia_changed=find(abs(IA_fin-(1:numel(IA_fin)))>0);
if(numel(  intersect(IA, dia_changed)  )>0)
    IA=IA_fin(IA)';
    fprintf(' [');fprintf('%3.0f',IA); fprintf(']');fprintf(' ==> ');
else
fprintf(' --> ');
end
IA=IA_fin(IA);
IB=DIAGN_star(IA);
IC=HDIAGN_star(IA);
fprintf('%10.0f ',C);
fprintf(' %s ',HDIAGN_star{IA});
fprintf('\n');

end
