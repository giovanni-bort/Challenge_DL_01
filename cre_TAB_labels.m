function [TAB,indici]=cre_TAB_labels(Labels,classes);

n_classes=numel(classes);
n_records=numel(Labels);fprintf('n_labels:%8.0f\n',n_records);
TAB=[];TAB(n_classes)=0;indici=[];
for ii=1:n_classes
    fprintf('cretab---> %6.0f %s --\n',ii,classes{ii});
      IIVal=find(Labels==classes{ii});
    indici(ii).ind=IIVal;
      TAB(ii)=TAB(ii)+numel(IIVal);

end

