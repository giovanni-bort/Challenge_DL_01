function [ind_TAB,ind_TAB_new,ind_new]=extend_train(Labels,classes)

fprintf(' Extend Training set - initail Labels:%6.0f\n',numel(Labels));
[TAB_ini,ind_TAB]=cre_TAB_labels(Labels,classes);
for iii=1:numel(TAB_ini), if(TAB_ini(iii)>0),fprintf('%3.0f TAB: %6.0f\n',iii,TAB_ini(iii)); end;end
ind_TAB_new=ind_TAB;
med_lab=median(TAB_ini(TAB_ini>0));
maxtab=sort(TAB_ini,'descend');
fprintf('Labels - median size(MS): %6.0f   2nd_label(2L)=  %6.0f  (=Total) (*** max 10 times)\n',med_lab,maxtab(2));
lab_new=[];
ind_new=[];lab_new=[];
n_TOT=0;
fprintf(' Label ->      N     MS/N     2L/N K1   K2       K1*N  K2*N = (K1+K2)N ->  Total       Ind1    Ind2   Ind1+2\n');
for ii=find(TAB_ini>0)
    N0=TAB_ini(ii);
    fprintf('%6.0f -> %6.0f',ii,TAB_ini(ii));
    if(TAB_ini(ii)<med_lab)
        KK=med_lab/TAB_ini(ii);fprintf(' * %6.2f',KK);
    else
        fprintf('        .');
    end
    if(TAB_ini(ii)<maxtab(2))
        KK2=(maxtab(2)/TAB_ini(ii))-1;fprintf(' * %6.2f',KK2);
    else
        KK2=0;
        fprintf('        .');
    end
    N1T=0;
    N2T=0;
    if(KK2>0),
        N1=fix(KK2);
        N2=KK2-N1;
        N1=min(N1,10); % *** MODIFICA *** max 10 times ****
        fprintf('%3.0f%6.2f',N1,N2);
        N1T=N0*N1;
        N2T=fix(N0*N2);
        fprintf('-> %6.0f%6.0f  = %6.0f   -> %6.0f',N1T,N2T,N1T+N2T,N0+N1T+N2T);
        new_ind=ind_TAB(ii).ind;
        for jj=1:N1
            new_ind=[new_ind  ;   ind_TAB(ii).ind];
        end
        ind22=ind_TAB(ii).ind(1:N2T);
        fprintf(' size%4.0f%2.0f ',size(ind_TAB(ii).ind(1:N2T)));
         fprintf('%5.0f%2.0f ',size(new_ind));
       new_ind=[new_ind  ;     ind_TAB(ii).ind(1:N2T)];
         fprintf('%5.0f%2.0f ',size(new_ind));
       
        ind_TAB_new(ii).ind=new_ind;
        
    end
    n_TOT=n_TOT+N0+N1T+N2T;
    fprintf('\n');
    
    ind_new=[ind_new ; ind_TAB_new(ii).ind];
    lab_new=[lab_new ];
    
end
    fprintf('Final total labels:%6.0f\n',n_TOT);
end

