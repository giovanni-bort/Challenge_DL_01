function [TAB_val,TAB_trn]=check_accuracy(Validation,Output,Classes,Training,Trained);

TAB_trn=0;
fprintf('Classes  ');
if(nargin>4),fprintf('( ------Training set------)');    end
fprintf('    ---------Validation set-------\n');
CLS=[];CLS(9,3)=0;
for ii=1:numel(Classes)
    fprintf(' %6s ',Classes{ii});
    if(nargin>4)
        IITing=find(Training==Classes{ii});
        IITrnd=find(Trained==Classes{ii});
        IITok=find(Trained(IITing)==Classes{ii});
        fprintf(' (%4.0f:%6.0f /%4.0f %6.1f%%)',numel(IITrnd),numel(IITok),numel(IITing),100*numel(IITok)/numel(IITing));
           CLS(ii,1)=100*numel(IITok)/numel(IITing);
           TAB_trn(ii,1)=numel(IITrnd);
           TAB_trn(ii,2)=numel(IITok);
           TAB_trn(ii,3)=numel(IITing);
 
    end
    IIVal=find(Validation==Classes{ii});
    IIP=find(Output==Classes{ii});
    IIPok=find(Output(IIVal)==Classes{ii});
    fprintf('%10.0f:%6.0f /%4.0f %6.1f%%',numel(IIP),numel(IIPok),numel(IIVal),100*numel(IIPok)/numel(IIVal));
    CLS(ii,2)=100*numel(IIPok)/numel(IIVal);
    fprintf('\n');
    TAB_val(ii,1)=numel(IIP);
    TAB_val(ii,2)=numel(IIPok);
    TAB_val(ii,3)=numel(IIVal);
    
end

% fprintf('Mean values: '); fprintf('%22.2f%%%31.2f%%',mean(CLS(:,1)),mean(CLS(:,2)));   fprintf('\n');
sum_trn=0;mean_trn=0;
sum_val=sum(TAB_val(:,3));

mean_val=mean(TAB_val(:,2)./TAB_val(:,3))*100;

fprintf('Mean/tot values: ')
if(nargin>4), sum_trn=sum(TAB_trn(:,3));mean_trn=mean(TAB_trn(:,2)./TAB_trn(:,3))*100;
    fprintf('%10.0f%7.1f%%          ',sum_trn,mean_trn);
end

    
fprintf(' %13.0f%7.1f%%\n',sum_val,mean_val);
