% prova-star



DIAGN_star=[270492004 	164889003 	164890007 	426627000 	713427006 	713426002 	445118002   39732003 	164909002 ...
            251146004 	698252002 	10370003 	284470004 	427172004  164947007 	111975006 	164917005 	47665007  ...
            59118001 	427393009 	426177001   426783006 	427084000 	63593006 	164934002 	59931005 	17338001 9999999 ];
HDIAGN_star ={ 'IAVB' ;   'AF'    ; 'AFL'     ; 'Barady'  ; 'CRBBB'    ; 'IRBBB'  ; 'LAnFB'   ; 'LAD'     ; 'LBBB' ; ...
               'LQRSV' ; 'NSIVCB' ; 'PR'      ;  'PAC'    ;  'PVC'     ;  'LPR'    ;  'LQT'   ; 'QaB'     ;  'RAD' ; ...
               'RBBB'  ;  'SA'    ; 'SB'      ; 'SNR'     ; 'STach'    ; 'SVPB'   ; 'TAb'     ; 'Tinv'    ; 'VEB' ; '*Other*' } ;

           
           
           
           
           fprintf('-------------  LABELS star DEFINITE  -------------------------------------\n');
        for i=1:numel(DIAGN_star)
             fprintf('%3.0f  ',i);
             fprintf('  %s ' ,HDIAGN_star{(i)});
             fprintf('  %15.0f ' ,DIAGN_star((i)));
             fprintf('\n');
         end
                   fprintf('-------------  -LABELS presenti- -------------------------------------\n');
                 
         DIAGN_STAR=[];for i=1:27,DIAGN_STAR{i}=num2str(DIAGN_star(i));end
         
         [I1,I2,I3]=intersect(classes,DIAGN_STAR);
         
         for i=1:numel(I1)
             fprintf('%3.0f %20s ',i,I1{i});
             fprintf('%6.0f',I2(i));             
             fprintf('%6.0f',I3(i));
             fprintf('  %s ' ,HDIAGN_star{I3(i)});
             fprintf('\n');
         end
         
         UNI_class=I1;
         fprintf('------------- MISSING LABELS -------------------------------------\n');
         [I1,I2,I3]=setxor(DIAGN_STAR,UNI_class);
         [I2S,IIS]=sort(I2);
         for i=1:numel(I1)
             fprintf('%3.0f %20s ',i,I1{i});
             fprintf('%6.0f',I2(i));         
             fprintf('  %s ' ,HDIAGN_star{I2(i)});
             fprintf('  %15.0f ' ,DIAGN_star(I2(i)));
             
%              fprintf('%6.0f',I3(i));
             fprintf('\n');
         end
 
                  fprintf('------------- MISSING LABELS -------------------------------------\n');
         [I1,I2,I3]=setxor(DIAGN_STAR,UNI_class);
         [I2S,IIS]=sort(I2);
         for i=1:numel(I1)
             fprintf('%3.0f %20s ',i,I1{IIS(i)});
             fprintf('%6.0f',I2(IIS(i)));         
             fprintf('  %s ' ,HDIAGN_star{I2(IIS(i))});
             fprintf('  %15.0f ' ,DIAGN_star(I2(IIS(i))));
             
%              fprintf('%6.0f',I3(i));
             fprintf('\n');
         end
 