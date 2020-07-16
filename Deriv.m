function[Zero1]=Deriv(Lead,Iz,J)
l=length(Lead);
Der(1:l-2)= Lead(1:l-2)-Lead(3:l);
   
% Zero crossing
l=length(Der);
[M Zero]=find(Der(1:l-1)<0&Der(2:l)>=0 | Der(1:l-1)>0&Der(2:l)<=0);
Zero=Zero+1;  % x correction
j=1;
Zero1=1;
for i=1:length(Zero)
    if Zero(i)>Iz & Zero(i)<J
        Zero1(j)=Zero(i);
        j=j+1;
    end
  end
    