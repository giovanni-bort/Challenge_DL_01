% Cross correlation
function[FP_,matrix]=cross_corr(lead,FP,Hz,Left,Right);
% function corrected 15.04.20
if FP(1)==FP(2); FP=FP(2:length(FP)); end
On=FP-round(Left*Hz/1000);    %80 ms left
Off=FP+round(Right*Hz/1000);  %100 ms right
if On(1)>0; k=1; else k=2; end
if Off(length(Off))>length(lead); l=2; else l=1; end

for i=k:length(FP)-l;
    s1=lead(On(i):Off(i));
    s2=lead(On(i+1):Off(i+1));
    X1=xcorr(s1,s2); %compute cross-correlation between vectors s1 and s2
    [m,FP_(i)]=max(X1);
    FP_(i)=FP(i)+FP_(i)-length(s1); %shift index FP_, as length(X1)=2*N-1; where N is the length of the signals
end
if k==1; FP_=[FP(1) FP_]; %to equalise the length
else FP_=[FP(1) FP_(2:length(FP_))];
end

% Matrix concatenation
Size=size(lead);        % raw or coloumn ?
if Size(1)~=1; lead=lead'; end
if FP_(1)>FP_(2); m=FP_(1); FP_(1)=FP_(2); FP_(2)=m;  end

% Boundery effects
m=2;matrix=[];
% % % fprintf('Left=%6.0f ',round(Left*Hz/1000));
% % % fprintf('lead:%6.0f FP_(%3.0f):',numel(lead),numel(FP_));fprintf('%5.0f',FP_):fprintf('\n');
if FP_(1)-round(Left*Hz/1000)>1; matrix=[lead(FP_(1)-round(Left*Hz/1000) : FP_(1)+round(Right*Hz/1000))];
else
    if FP_(2)-round(Left*Hz/1000)>1; matrix=lead(FP_(2)-round(Left*Hz/1000) : FP_(2)+round(Right*Hz/1000));
    else
        if(FP_(3)-round(Left*Hz/1000)>1), matrix=lead(FP_(3)-round(Left*Hz/1000) : FP_(3)+round(Right*Hz/1000)); m=3;end
    end
end
if FP_(length(FP_))+Right>length(lead); FP_=FP_(1:length(FP_)-1); end
for i=m:length(FP_);
  if(FP_(i)-round(Left*Hz/1000)<1),
      fprintf('***i:%6.0f -> %6.0f%6.0f ',i,FP_(i),round(Left*Hz/1000));
      fprintf('  matrix: %6.0f%6.0f  %6.0f%6.0f\n',size(matrix),round(Left*Hz/1000),round(Right*Hz/1000));
      new(1:round(Left*Hz/1000)+round(Right*Hz/1000)+1)=0;
  else %********** MODIFY *******
    new=lead(FP_(i)-round(Left*Hz/1000) : FP_(i)+round(Right*Hz/1000));
  end
   matrix=[matrix; new];
end
%fprintf(' **** cross_corr: matrix:%6.0f%6.0f\n',size(matrix));