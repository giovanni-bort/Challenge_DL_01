
%================================================================================================================
%================================== INI program of classification ===============================================    

if(Freq==500),fprintf('frequency:500Hz\n');
else
    ECG=Interpolation(ECG',Freq,500);
    ECG=ECG';
   fprintf('freq:%6.0f  -> resampled at 500Hz   size ECG:%4.0f%8.0f\n', Freq,size(ECG))

end
   Hz=500;



    ECG=ECG/1000;
    I=ECG(1,:); II=ECG(2,:); III=ECG(3,:);
    aVR=ECG(4,:); aVL=ECG(5,:); aVF=ECG(6,:);
    V1=ECG(7,:); V2=ECG(8,:); V3=ECG(9,:);
    V4=ECG(10,:); V5=ECG(11,:); V6=ECG(12,:);

    
    
   [A,a]=find(I>2); I(a)=2; [A,a]=find(II>2); II(a)=2;
   [A,a]=find(III>2); III(a)=2;[A,a]=find(V1>2); V1(a)=2;
   [A,a]=find(V2>2); V2(a)=2;[A,a]=find(V3>2); V3(a)=2;
   [A,a]=find(V4>2); V4(a)=2; [A,a]=find(V5>2); V5(a)=2;
   [A,a]=find(V6>2); V6(a)=2;

    [A,a]=find(I<-2); I(a)=-2; [A,a]=find(II<-2); II(a)=-2;
    [A,a]=find(III<-2); III(a)=-2;[A,a]=find(V1<-2); V1(a)=-2;
    [A,a]=find(V2<-2); V2(a)=-2;[A,a]=find(V3<-2); V3(a)=-2;
    [A,a]=find(V4<-2); V4(a)=-2; [A,a]=find(V5<-2); V5(a)=-2;
    [A,a]=find(V6<-2); V6(a)=-2;

PVC(j)=0;
OldMI(j)=0; % old miocardial infarction
AnMI(j)=0;  % anterior miocardial infarction
MI(j)=0;    % miocardial infarction
LAnFB(j)=0;
SA(j)=0; 
SNR(j)=0; 

% QRS detectin
    d1=I(1:2:length(II));      % 250 Hz
    d2=V2(1:2:length(II));
        if(sum(abs(d1))==0), d1=II (1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d1))==0), d1=III(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d1))==0), d1=aVR(1:2:length(II)); end        %MODFIED 22.04.20
        if(sum(abs(d2))==0), d2=V3(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d2))==0), d2=V4(1:2:length(II)); end       %MODFIED 22.04.20
        if(sum(abs(d2))==0), d2=V5(1:2:length(II)); end       %MODFIED 22.04.20
    
    [QRS]=QRS_det(d1,d2);      % go to function
    QRS=QRS*2;  % 500Hz
     figure(1); clf
      x=(1:1:length(I))/Hz;
      plot(x,V1, QRS/Hz, V1(QRS), 'ro')

if QRS(1)==QRS(2); QRS=QRS(2:length(QRS)); end
if length(I)/length(QRS)>=1000
ecg=II;
if(sum(abs(ecg(:)))==0),ecg=III;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V1;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V2;end    %MODFIED 22.04.20
if(sum(abs(ecg(:)))==0),ecg=V3;end    %MODFIED 22.04.20


fs=1000;
QRS=0;
[QRS,sign,en_thres] = qrs_detect2(ecg',0.25,0.6,fs);
figure(1); clf
x=(1:1:length(I))/Hz;
plot(x,II, QRS/Hz, II(QRS), 'ro')
end

I_AVB(j)=0; AFL(j)=0; AF(j)=0; PVC(j)=PVC(j); PAC(j)=0; LBBB(j)=0;
CRBBB(j)=0; IRBBB(j)=0;Normal(j)=0; STD(j)=0; STE(j)=0;

if mean(I)~=0
%======================================================================
%		Preprocessing
%======================================================================n=7;					% 7*2 = 14 ms smoothing interval on each side
T=1/Hz;					% [s] sampling period
Fc=0.64;				% [Hz]
a=1; b=[1 1 1 1 1 1 1 1 1 1]/10;
n=7;

 X=I;        Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	I=Y;	% Drift filtration
 
 X=II;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	II=Y;	% Drift filtration
 
 X=V1;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V1=Y;	% Drift filtration
 
 X=V2;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V2=Y;	% Drift filtration
 
 X=V3;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V3=Y;	% Drift filtration
 
 X=V4;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V4=Y;	% Drift filtration
 
 X=V5;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V5=Y;	% Drift filtration
 
 X=V6;       Y=filter(b,a,X);        % 50 Hz
 X=Y;        [Y]=smooth(X,n);			% Smoothing
 X=Y;		[Y]=drift_Ivo(X,Fc,T);	V6=Y;	% Drift filtration
 
% Ventricular fibrillation ?
VF(j)=0;
Lead=I; Iz=1; J=length(Lead);
    [Zero1]=Deriv(Lead,Iz,J);
if length(Zero1)/length(Lead)<0.019
    VF(j)=1;
end

if isempty(QRS)==0
% Fiducial point = at the max peak
FP=0;
left=80*Hz/1000;    %180 ms
right=100*Hz/1000;   %180 ms

if QRS(1)-left<=0; QRS=QRS(2:length(QRS)); end
if QRS(1)-left<=0; QRS=QRS(2:length(QRS)); end

Max=max(II(QRS(1)-left:QRS(1)+right)); Min=abs(min(II(QRS(1)-left:QRS(1)+right)));
if QRS(1)- left <= 0; QRS=QRS(2:length(QRS)); end
if Max>=Min;
    for i=1:length(QRS)-2;
        [M,m]=max(II(QRS(i)-left:QRS(i)+right));
        FP(i)=QRS(i)-left+m;
    end    
else
    for i=1:length(QRS)-2;
        [M,m]=min(II(QRS(i)-left:QRS(i)+right));
        FP(i)=QRS(i)-left+m ;
    end    
end


% PAC(j)
PAC(j)=0; 
if FP(1)==FP(2); FP=FP(2:length(FP)); end
RR=FP(2:length(FP))-FP(1:length(FP)-1);
if std(RR)>35     PAC(j)=1; end


% PVC
Lead=II;

if(sum(abs(Lead(:)))==0),Lead=II;end   %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=III;end  %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=V1;end   %*** MODIFIED 22.4.20
if(sum(abs(Lead(:)))==0),Lead=V2;end   %*** MODIFIED 22.4.20

plot(x,Lead, FP/Hz,Lead(FP), 'ro');
if mean(Lead(FP))>0
    if max(Lead(FP))*0.50>mean(Lead(FP))
    PVC(j)=1; %PAC(j)=0;
    end
    if min(Lead(FP))<0.50*mean(Lead(FP))
    PVC(j)=1; %PAC(j)=0;
    end
else
    if min(Lead(FP))*0.50<mean(Lead(FP))
    PVC(j)=1; %PAC(j)=0;
    end
    if max(Lead(FP))>0.50*mean(Lead(FP))
    PVC(j)=1; %PAC(j)=0;
    end
end

%STach sinus tachycardia
RRm=mean(RR);
if PAC(j)==0 && PVC(j)==0 && 1000/Hz*RRm<600;
    STach(j)=1;
else STach(j)=0;
end

% Brady
if PAC(j)==0 && PVC(j)==0 && 1000/Hz*RRm>1000;
    Brady(j)=1;
else Brady(j)=0;
end

% Sinus arrhythmia
if PAC(j)==1 && PVC(j)==0;
    SA(i)=1;
else SA(i)=0;
end

% Sinus rhythm
if PAC(j)==0 && PVC(j)==0
    SNR(j)=1;
else SNR(j)=0;
end


% SB
if Brady(j)==1 && SNR(j)==1; SB(j)=1; else SB(j)=0; end
    

RRm=mean(RR);
Left=650*Hz/1000;    %650 ms left
Right=round(RRm)+150*Hz/1000;
Right=round(Right/2)*2; % Right must be even

if RRm<250
    Left=600*Hz/1000;
    Right=Left+round(RRm)+50*Hz/1000;
    Right=round(Right/2)*2;
end


% Average of all leads
A=find(RR<150);
FP(A)=0;

for i=length(A):-1:1;
    n=A(i);
    FP=[FP(1:n-1), FP(n+1:length(FP))];
end
if FP(length(FP))==0; FP=FP((1):length(FP)-1); end

% exclude FP of PVC
if PAC==1;
    k=1;
    B=mean(II(FP));
    for i=1:length(FP)
        if abs(B-II(FP(i)))<0.5*abs(B)
            A(k)=FP(i); k=k+1;
        end
    end
FP=A;
end
        

figure(1)
clf
plot(x,II, FP/Hz,II(FP), 'ro')
title(num2str(j))


[FP_,matrix]=cross_corr(I,FP,Hz,Left,Right);     I_matrix=matrix;  I_mean=mean(I_matrix);
[FP_,matrix]=cross_corr(II,FP,Hz,Left,Right);    II_matrix=matrix; II_mean=mean(II_matrix);
[FP_,matrix]=cross_corr(III,FP,Hz,Left,Right);   III_matrix=matrix; III_mean=mean(III_matrix);
if mean(V1)~=0      %Peripheral leads only
[FP_,matrix]=cross_corr(aVF,FP,Hz,Left,Right);   aVF_matrix=matrix; aVF_mean=mean(aVF_matrix);
[FP_,matrix]=cross_corr(aVR,FP,Hz,Left,Right);   aVR_matrix=matrix; aVR_mean=mean(aVR_matrix);
[FP_,matrix]=cross_corr(V1,FP,Hz,Left,Right);    V1_matrix=matrix; V1_mean=mean(V1_matrix);
[FP_,matrix]=cross_corr(V2,FP,Hz,Left,Right);    V2_matrix=matrix; V2_mean=mean(V2_matrix);
[FP_,matrix]=cross_corr(V3,FP,Hz,Left,Right);    V3_matrix=matrix; V3_mean=mean(V3_matrix);
[FP_,matrix]=cross_corr(V4,FP,Hz,Left,Right);    V4_matrix=matrix; V4_mean=mean(V4_matrix);
[FP_,matrix]=cross_corr(V5,FP,Hz,Left,Right);    V5_matrix=matrix; V5_mean=mean(V5_matrix);
if mean(V6)~=0
[FP_,matrix]=cross_corr(V6,FP,Hz,Left,Right);    V6_matrix=matrix; V6_mean=mean(V6_matrix);
end
end
x1=1/Hz:1/Hz:length(I_mean)/Hz;
figure(2); clf
subplot(3,3,1); plot(x1,I_mean); title('I')
subplot(3,3,2); plot(x1,II_mean); title('II')
subplot(3,3,3); plot(x1,III_mean); title('III')
if mean(V1)==0; 
    V1_mean(1:length(I_mean))=0;
    V2_mean(1:length(I_mean))=0;
    V3_mean(1:length(I_mean))=0;
    V4_mean(1:length(I_mean))=0;
    V5_mean(1:length(I_mean))=0;
    V6_mean(1:length(I_mean))=0;
end
    if mean(V6)==0;     V6_mean(1:length(I_mean))=0; end

x2=(1:1:length(V1_mean))/Hz;
subplot(3,3,4); plot(x2,V1_mean); title('V1')
subplot(3,3,5); plot(x2,V2_mean); title('V2')
subplot(3,3,6); plot(x2,V3_mean); title('V3')
subplot(3,3,7); plot(x2,V4_mean); title('V4')
subplot(3,3,8); plot(x2,V5_mean); title('V5')
x3=(1:1:length(V6_mean))/Hz;
subplot(3,3,9); plot(x3,V6_mean); title('V6')

% Find the izoelectric point in II or V1
a=90;  b=200;  %  Offset
a=fix(2*90*Hz/1000);  b=fix(2*200*Hz/1000);  %  Offset
Flat=20*Hz/1000; 
Lead=I_mean;
b=min(b,numel(Lead));
[Max,t1]=max(Lead(a:b));
[Min,t2]=min(Lead(a:b));
Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end

To=max(From-60, 12);
[Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
Iz_1=Pnt-Flat;

Flat=20*Hz/1000;        % 20 ms
From=From+30;            			% from the point of the peak
To=From+80*Hz/1000;                % to 100 ms right ...
[Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
J_1=Pnt;

x1=(1:1:length(I_mean))/Hz;
subplot (3,3,1);
plot(x1,Lead), hold on
plot(Iz_1/Hz, Lead(Iz_1), 'ro', J_1/Hz, Lead(J_1), 'ro')
title ('I')

Lead=V1_mean;
Iz2=Iz_1;

if V1_mean~=0;
[Max,t1]=max(Lead(a:b));
[Min,t2]=min(Lead(a:b));
Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end
From1=From;
To=max(From-60, 12);
[Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function

Iz_2=Pnt-Flat;       % 20 ms
From=From+30;            			% from the point of the peak
To=From+80*Hz/1000;                % to 100 ms right ...
[Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
J_2=Pnt;

x1=(1:1:length(Lead))/Hz;
subplot (3,3,4);
plot(x1,Lead, Iz_2/Hz, Lead(Iz_2), 'ro', J_2/Hz, Lead(J_2), 'ro')
    Iz=max(Iz_1, Iz_2);

Iz2=Iz;
if abs(J_1-J_2) <20; J=max(J_1, J_2);
    else J=round((J_1+J_2)/2);
end


% RBBB(j) criteria 
% 'M' shaped  QRS in V1; wide S in V6 & I
%	Searching for QRS segmentation 
RBBB1=0; RBBB2=0; RBBB3=0; RBBB4=0; CRBBB(j)=0; IRBBB(j)=0; LBBB(j)=0;
a=15;
    
Lead=V2_mean;
From=From1;
To=max(From-60, 12);
[Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
Iz=Pnt-Flat;

    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz
       
    
    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,5);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    plot(Zero1/Hz, Lead    (Zero1), 'go')   
    if length(Zero1)>=3 %& no==0; 
        title ('RBBB'); RBBB1=1; 
    else title ('V2'); RBBB1=0; end


    Lead=V1_mean;
    From=From1;
    To=max(From-60, 12);
    [Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
    Iz=Pnt-Flat;

    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz

    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,4);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko');
    hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3 %& no==0;
        title ('RBBB'); RBBB2=1; 
    else title ('V1'); RBBB2=0; end

     Lead=V3_mean;
    From=From1;
    To=max(From-60, 12);
    [Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
    Iz=Pnt-Flat;
    
    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function

    x1=(1:1:length(V1_mean))/Hz;
    Iz=Iz-a;   J=J+a;       % restore J & Iz
    subplot(3,3,6);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko',J/Hz,Lead (J), 'ko'); hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3 %& no==0; 
        title ('RBBB'); RBBB3=1; 
    else title ('V3'); RBBB3=0; end
    
    
    Lead=V4_mean;
    From=From1;
    To=max(From-60, 12);
    [Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
    Iz=Pnt-Flat;
    Iz=Iz+a;    J=J-a;
    [Zero1]=Deriv(Lead,Iz,J); % go to function
    Iz=Iz-a;   J=J+a;   % restore J & Iz
    
    x1=(1:1:length(V1_mean))/Hz;
    subplot(3,3,7);
    plot(x1,Lead, Iz/Hz,Lead (Iz), 'ko', J/Hz,Lead(J),'ko'); hold on
    plot(Zero1/Hz, Lead(Zero1), 'go')   
    if length(Zero1)>=3 %& no==0; 
        title ('RBBB'); RBBB4=1; 
    else title ('V4'); RBBB4=0; end
    
RBBB(j)=0;
if RBBB1|RBBB2|RBBB3|RBBB4~=0;
    RBBB(j)=1;
    if J-Iz<120*Hz/1000; IRBBB(j)=1; CRBBB(j)=0; end
    if J-Iz>120*Hz/1000; IRBBB(j)=0; CRBBB(j)=1; end
end
    


% LBBB(j) criteria 
% monophasic R, no Q in I & V6; QS in V1
LBBB(j)=0;

if J-Iz>120*Hz/1000    
    Lead=I_mean; Minus=0; Plus=0;
    Lead=Lead-Lead(Iz);
    for i=Iz:J
        if Lead(i)>0; Plus=Plus+Lead(i);
        else Minus=Minus+Lead(i);
        end
    end

    if Plus > 15*abs(Minus); LBBB(j)=1;
    %elseif Minus> 5*Plus; LBBB(j)=0;
    end
    
    Lead=V6_mean; 
    a=90;  b=200;  %  Offset
    a=fix(2*90*Hz/1000);  b=fix(2*200*Hz/1000);  %  Offset
    b=min(b,numel(Lead));
    [Max,t1]=max(Lead(a:b));
    [Min,t2]=min(Lead(a:b));
    Min=abs(Min);
if Max>=Min; 
    Crit=0.02*Max;	% make Crit proportional to the QRSampl;
    From=t1+a;
else Crit=0.02*Min;	% make Crit proportional to the QRSampl;j
    From=t2+a;
end

    prag=0.8*Lead(From);
    for i=From:-1:From-40
        if Lead(i)<=prag; break
        end
    end
    pl=i;
    
    for i=From:From+40
        if Lead(i)<=prag
            0.8*Lead(From); break
        end
    end
    pr=i;
    
    if pr-pl>70*Hz/1000; LBBB(j)=1;  end
    if pr-pl<30*Hz/1000;LBBB(j)=0;   end
end
  
   figure(2);
   Lead=II_mean;
   subplot(3,3,2);
    x1=(1:1:length(Lead))/Hz;
    plot(x1,II_mean, Iz/Hz, II_mean(Iz), 'ko',  J/Hz, II_mean(J), 'ko')
    title ('II')
   
    Lead=III_mean;
   subplot(3,3,3);
    x1=(1:1:length(Lead))/Hz;
    plot(x1,III_mean, Iz/Hz, Lead(Iz), 'ko',  J/Hz, Lead(J), 'ko')
    title ('III')
    
    Lead=V5_mean;
    subplot(3,3,8);
    plot(x1,Lead, Iz/Hz, Lead(Iz), 'ko',  J/Hz, Lead(J), 'ko')
    title ('V5')

if length(Lead)==length(V6_mean);
    Lead=V6_mean;
    subplot(3,3,9);
    plot(x1,Lead, Iz/Hz, Lead(Iz), 'ko',  J/Hz, Lead(J), 'ko')
    title ('V6')
end

Iz2=Iz; J2=J;

% STD(j) STE(j)
ST1=I_mean(J2)-I_mean(Iz2); 
ST2=II_mean(J2)-II_mean(Iz2);
ST3=III_mean(J2)-III_mean(Iz2);
ST4=V1_mean(J2)-V1_mean(Iz2);
ST5=V2_mean(J2)-V2_mean(Iz2);
ST6=V3_mean(J2)-V3_mean(Iz2);
ST7=V4_mean(J2)-V4_mean(Iz2);
ST8=V5_mean(J2)-V5_mean(Iz2);
ST9=V6_mean(J2)-V6_mean(Iz2);

ma=[ST1,ST2,ST3,ST4,ST5,ST6,ST7,ST8,ST9];
[Max,m]=max(abs(ma));
ST=ma(m);
if ST>0.1; STE(j)=1; else STE(j)=0; end
if ST<-0.1 STD(j)=1; else STD(j)=0; end

% AMIs - acute miosardial ischemia
if STE(j)==1; AMIs(j)=1; else AMIs(j)=0; end


% old miocardial infarction
% anterior miocardial infarction
% miocardial infarction
% AnMI
if m==1; Lead=I_mean;
elseif m==2; Lead=II_mean;
elseif m==3; Lead=III_mean;
elseif m==4; Lead=V1_mean;
elseif m==5; Lead=V2_mean;
elseif m==6; Lead=V3_mean;
elseif m==7; Lead=V4_mean;
elseif m==8; Lead=V5_mean;
elseif m==9; Lead=V6_mean;
end

if ST>0.2
   if Lead(J2:J2+70)>0
        AnMI(j)=1;
   end
end

if ST<-0.2
    if Lead(J2:J2+70)<0
        AnMI(j)=1;
    end
end


% MI
ma=[ST1,ST2,ST3,ST4,ST5,ST6,ST7,ST8,ST9];
[Max,m]=max(abs(ma));
ST=ma(m);
if m==1; Lead=I_mean;
elseif m==2; Lead=II_mean;
elseif m==3; Lead=III_mean;
elseif m==4; Lead=V1_mean;
elseif m==5; Lead=V2_mean;
elseif m==6; Lead=V3_mean;
elseif m==7; Lead=V4_mean;
elseif m==8; Lead=V5_mean;
elseif m==9; Lead=V6_mean;
end
    
if ST>0.1
   if Lead(J2:J2+70)>0
       MI(j)=1;
   end
 end

if ST<-0.1
   if Lead(J2:J2+70)<0
       MI(j)=1;
end
 
end


% Old miocardial infarction
%  Q wave in lead III wider than 1 mm and
%  Q wave in lead aVF wider than 0.5 mm and
%  Q wave of any size in lead II

 k1=0; k2=0; k3=0;
 Iz=Iz+20; J=J-10;
 Lead=III_mean;
     [Zero1]=Deriv(Lead,Iz,J);
%      figure(2)
%      subplot(3,3,3); hold on
%      plot(Zero1/Hz, Lead(Zero1), 'ro')
 if Lead(Zero1(1))<0
     for i=Iz:Iz+30
         if Lead(i)<0
             k1=k1+1;
         else break
         end
     end
 end
     
Lead=aVF_mean;
  [Zero1]=Deriv(Lead,Iz,J);
if Lead(Zero1(1))<0
    for i=Iz:Iz+30
        if Lead(i)<0
             k2=k2+1;
        else break
        end
     end
end


Lead=II_mean;
[Zero1]=Deriv(Lead,Iz,J);
%figure(2)
%subplot(3,3,2); hold on
%plot(Zero1/Hz, Lead(Zero1), 'ro')
if Lead(Zero1(1))<0
     for i=Iz:Iz+30
         if Lead(i)<0
             k3=k3+1;
         else break
         end
     end
 end
 
if k1>16 && k2>8 && k3>4
     OldMI(j)=1;
end

 
%---------
% I_AVB(j)
% Find P wave in Average ecg

% Positive P-wave


prag2=0.0050;     % 5 uV threshold

P_wave=0; 
Lead=II_mean;

for i=Iz2-20:-1:1;
    if Lead(i+20)-Lead(i)>4*prag2 && Lead(i+20)-Lead(i+40)>1*prag2 &&...
            Lead(i+20)-Lead(i+10)>prag2 && Lead(i+20)-Lead(i+30)>prag2/2;
        P_wave=1;
        break
    end
end

P_x=i+20;

if P_wave==0
    Lead=I_mean;
     [M,m]=max(Lead(20:Iz2));
     f1=find(M-Lead(20:m)>prag2*5);
     f2=find(M-Lead(m:Iz2-20)>prag2*5);
 if isempty(f1)==0 && isempty(f2)==0
     P_wave = 1;
     P_x=m+20;
end
end


if P_wave==0            % negative P-wave
     [M,m]=min(Lead(20:Iz2));
     f1=find(Lead(20:m)-M>prag2*5);
     f2=find(Lead(m:Iz2-20)-M>prag2*5);
 if isempty(f1)==0 && isempty(f2)==0
     P_wave = 1;
     P_x=m+20;
 end

end

Iz=Iz2;

if P_x<30; P_wave=0; end;% too far
if Iz-P_x<40; P_wave=0; end % too close
    
Lead=II_mean;

figure(2);
    x=(1:1:length(I_mean))/Hz;
    subplot(3,3,2);
    
    
    if P_wave==1
    plot(x, Lead, P_x/Hz, Lead(P_x), 'go', Iz/Hz, Lead(Iz), 'ko', J/Hz, Lead(J), 'ko');
    else plot(x, Lead, Iz/Hz, Lead(Iz), 'ko',  J/Hz, Lead(J), 'ko');
    end
    title('II')    
    
if P_wave==1 && Iz-P_x > 125*Hz/1000; I_AVB(j)=1; LPR(j)=1;
else I_AVB(j)=0; LPR(j)=0;
end

Izsave=Iz;
Jsave=J;

% Double check for AF(j)
Lead=I;
AF(j)=0;
B=0;
for i=1:length(FP)-1
    RR(i)=FP(i+1) - FP(i);
end
[R,m]=max(RR(1:length(RR)-2));
if isempty(m)==0
    From=FP(m+1)-20;
    To=From-4;
    [Pnt,Crit]=Flat_l(Lead,Crit,From,To,Flat);  % go to function
    Iz_=Pnt-Flat;
    From=FP(m);
    From=From+20;            			% from the point of the peak
    To=From+100*Hz/1000;                % to 100 ms right ...
    [Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);  % go to function
    J_=Pnt;
   
    Iz=J_; J=Iz_;
    
     figure(1); hold on
     
  %      plot(J/Hz, Lead(J), 'ro', Iz/Hz, Lead(Iz), 'ro')

   LAnFB(j)=0;
  [Zero1]=Deriv(Lead,Iz,J); % go to function
    B=[Iz,Zero1,J]; 
           
    if P_wave==0 && PAC(j)==1;
        if length(B)>25; AF(j)=1; 
        elseif length(B)>19; AFL(j)=1;
        else LAnFB(j)=1;
        end
        figure(1); hold on
        plot(B/Hz, Lead(B), 'go')
    end
end

    % II-AVB
II_AVB(j)=0;
if I_AVB(j) ==1 && PVC(j)==0 && PAC(j)==1 
    [Ma,ma]=max(RR);
    if ma~=1; 
        if RR(ma)>1.9*RR(ma-1);
             II_AVB(j)=1;
        end
    end
end
         
end


%---------------------------
% Tend detection
% Searching for big peaks by 'wings' function
clear Wing;
s=30*Hz/1000;   % 40 ms    
Lead=V2_mean;
J=Jsave;

Wing(J+s:length(Lead)-s)=...
    (Lead(J:length(Lead)-2*s)-Lead(J+s:length(Lead)-s)).*...
    (Lead(J+s:length(Lead)-s)-Lead(J+2*s:length(Lead)));
Wing=Wing(J+s:length(Wing));

To=200;

[Mi2,mi2]=max(Wing(50:min([To*Hz/1000 length(Wing)])));
Tp2=J+mi2+s-1+50;

Lead=Lead-Lead(J);

[Mi1,mi1]=min(Wing(50:min([To*Hz/1000 length(Wing)])));
Tp1=J+mi1+s-1+50;

Tp=min([Tp1 Tp2]);

[Ma,ma]=max(abs(Lead(J+100*Hz/1000:length(Lead)-150*Hz/1000)));

Tp=J+100*Hz/1000+ma-1;

if Lead(Tp)>=0
    [Mi,mi]= min(Lead(Tp:length(Lead)-100*Hz/1000));
else
    [Mi,mi]= max(Lead(Tp:length(Lead)-100*Hz/1000));
end
T_end=Tp+mi-1;


% Second check
[Ma,ma]=max(Lead(J+100*Hz/1000:length(Lead)-50*Hz/1000));
[Mi,mi]=min(Lead(J+100*Hz/1000:length(Lead)-50*Hz/1000));
Tp1=Ma-Lead(J)+Ma-Lead(T_end);
Tp2=abs(Mi-Lead(J))+abs(Mi-Lead(T_end));
if Tp1>Tp2; Tp=ma+J+100*Hz/1000-1;
else Tp=mi+J+100*Hz/1000-1;
end

if Lead(Tp)>=0
    [Mi,mi]= min(Lead(Tp:length(Lead)-50*Hz/1000));
else
    [Mi,mi]= max(Lead(Tp:length(Lead)-50*Hz/1000));
end
T_end=Tp+mi-1;
figure(2)
subplot(3,3,5)
hold on
plot(T_end/Hz, V2_mean(T_end), 'go', Tp/Hz, V2_mean(Tp), 'go')

subplot(3,3,8)
plot(x1, aVR_mean, T_end/Hz, aVR_mean(T_end), 'go', Tp/Hz, aVR_mean(Tp), 'go',...
    Izsave/Hz, aVR_mean(Izsave), 'ko', Jsave/Hz, aVR_mean(Jsave), 'ko')
title('aVR')

% LQT
LQT(j)=Tp-Izsave;
if LQT(j)>340*Hz/1000; LQT(j)=1; else LQT(j)=0; end  

% T inversion
if V1_mean(Tp)<V1_mean(T_end) && ...
        V2_mean(Tp)<V2_mean(T_end) && ...
        V3_mean(Tp)<V3_mean(T_end); 
        TInv(j)=1;
else
    TInv(j)=0;
end


 % Search fot inflection point
c=length(Lead)-round(50*Hz/1000);
Tfrom=Tp+round(80*Hz/1000);

a=Lead(Tfrom+6:c-6)-Lead(Tfrom+3:c-9);
b=Lead(Tfrom+6:c-6)-Lead(Tfrom+9:c-3);
infl=a.*b;
inf=find(infl>=0);

if isempty(inf)
    T_end=T_end; 
else
    infl=inf(1);
    T_end=Tfrom+infl+6;
end


if AF(j)==1; I_AVB(j)=0; end

if PAC(j)~0; PAC(j)=1; end

% ST interal abnormality
if I_AVB(j) + AF(j) + PVC(j) + PAC(j) + LBBB(j) + CRBBB(j) ==0;
    STIAb(j) = 0; NSSTTA(j)=0;
else STIAb(j)=0; NSSTTA(j)=0;
end

% Conflictig are AF(j)-AVB, AF(j)-PAC(j), AF(j)-N, LBBB(j)-RBBB(j), AVB-RBBB(j) 
if AF(j)==1; I_AVB(j)=0; end
if LBBB(j)==1; CRBBB(j)=0; end 
if CRBBB(j)==1; LBBB(j)=0; end
%if OldMI(j)==1; STIAb(j)=0; NSSTTA(j)=0; end
%if VF(j)==1;
%     AF(j)=0;IAVB(j)=0; LBBB(j)=0; Nor(j)=0; PAC(j)=0; PVC(j)=0; 
%     CRBBB(j)=0; STD(j)=0; STE(j)=0; STIAb(j)=0; Brady(j)=0; IIAVB(j)=0; II_AVB(j)=0,...
%     AMIs(j)=0; STach(j)=0; OldMI(j)=0; SNR(j)=1; 
% else SNR(j)=0;
% end

% VCG
% QRS vector loop
%       Cut the records to get just the QRS
Iz=Izsave; J=Jsave;
Im=I_mean(Iz:J); IIm=II_mean(Iz:J); V1m=V1_mean(Iz:J); V2m=V2_mean(Iz:J);
V3m=V3_mean(Iz:J); V4m=V4_mean(Iz:J); V5m=V5_mean(Iz:J); V6m=V6_mean(Iz:J);

%-------------------------------------------------------------------------
%		Orthogonal lead evaluation
%       Dower’s transform:
% [X Y Z] = [V1 V2 V3 V4 V5 V6 I II]
%   X=[-0.172V1 -0.073V2 0.122V3 0.231V4 0.239V5 0.193V6 0.156I -0.009II]
%   Y=[0.057V1 -0.019V2 -0.106V3 -0.022V4 0.040V5 0.048V6 -0.227I 0.886II]
%   Z=[-0.228V1 -0.310V2 -0.245V3 -0.063V4 0.054V5 0.108V6 0.021I 0.102II]
x = -0.172*V1m -0.073*V2m +0.122*V3m +0.231*V4m +0.239*V5m +0.193*V6m +0.156*Im -0.009*IIm;
y = 0.057*V1m -0.019*V2m -0.106*V3m -0.022*V4m +0.040*V5m +0.048*V6m -0.227*Im +0.886*IIm; 				
z = -0.228*V1m -0.310*V2m -0.245*V3m -0.063*V4m +0.054*V5m +0.108*V6m +0.021*Im +0.102*IIm;
x = x-x(1); y=y-y(1); z=z-z(1);		% orth. leads - isoelectric point

% Close the loops
x=[x,0]; y=[y,0]; z=[z,0];

% Maximal vector
[QRS_vect_, max_pnt]= max(sqrt(x.^2+y.^2));
QRS_vect_=(round(QRS_vect_*100))/100;
QRSx_=x(max_pnt); QRSy_=-y(max_pnt);

% Angle of the QRS Vector
[QRS_ang_]=quadr(QRSx_,QRSy_);
QRS_ang_=round(QRS_ang_);

figure(3);
clf; subplot(1,2,1);
plot(x,y,'b','LineWidth',2)
xlabel('X ->'); ylabel('<- Y');
clf; subplot(1,2,1);
plot(x,y,'b','LineWidth',2)
xlabel('X ->'); ylabel('<- Y');
d_x=max(x)-min(x); d_y=max(y)-min(y); 
if d_x>=d_y; 
    x_left=min(x)-0.01; x_right=max(x)+0.01;
    y_down=min(y)-(d_x-d_y)/2+0.01; y_up=max(y)+(d_x-d_y)/2+0.01;
else
    y_down=min(y)-0.01; y_up=max(y)+0.01;
    x_left=min(x)-(d_y-d_x)/2+0.01; x_right=max(x)+(d_y-d_x)/2+0.01;
end


axis([x_left x_right y_down y_up]);
axis square;
view(0, -90);									% Y axis down

hold on
plot(0,0,'ro', QRSx_,-QRSy_,'ro')
plot([0,QRSx_],[0,-QRSy_], 'r')

RAD(j)=0; LAD(j)=0;
if QRS_ang_>=60 && QRS_ang_<=180; RAD(j)=0; LAD(j)=0; end
if QRS_ang_>=0 && QRS_ang_<=60; RAD(j)=0; LAD(j)=1; end
if QRS_ang_>=180 && QRS_ang_<=315; RAD(j)=1; LAD(j)=0; end
if QRS_ang_>=316 && QRS_ang_<=360; RAD(j)=0; LAD(j)=1; end

% PR
if P_wave==0 && PAC(j)==0 && (RAD(j)==1 || LAD(j)==1)
    PR(j)=1; AF(j)=0;
else PR(j)=0;
end

% QAb
QAb(j)=0;

% NSIVCB
if RBBB(j)==1 || LBBB(j)==1 || LQT(j) ==1; 
    NSIVCB(j)=1;
else NSIVCB(j)=0;
end

% LQRSV
Iz=Izsave; J=Jsave; LQRSV(j)=0;
[Ma1,ma1]=max(abs(I_mean(Iz:J)));
[Ma2,ma2]=max(abs(II_mean(Iz:J)));
[Ma3,ma3]=max(abs(III_mean(Iz:J)));
Mi=min([Ma1 Ma2 Ma3]);
if Mi<0.25; LQRSV(j)=1; end

% TAb
TAb(j)=0;


TEXT1=['AF=', num2str(AF(j)), '  AFL=', num2str(AFL(j)), '  VF=', num2str(VF(j))',...
    '  LAnFB=', num2str(LAnFB(j)), '  IAVB=', num2str(I_AVB(j))', '  IIAVB=', num2str(II_AVB(j)),...
    '  LBBB=', num2str(LBBB(j)), '  CRBBB=', num2str(CRBBB(j)), '  IRBBB=', num2str(IRBBB (j)),...
    '  PAC=', num2str(PAC(j)), '  PVC=', num2str(PVC(j)),...
    '  STD=', num2str(STD(j)), '  STE=', num2str(STE(j))];
    
TEXT2=['STIAb=', num2str(STIAb(j)), '  Brady=', num2str(Brady(j)), '  AMIs=', num2str(AMIs(j)),...
    '  NSSTTA=', num2str(NSSTTA(j)), '  STach=', num2str(STach(j)), '  SNR=', num2str(SNR(j)),...
    '  PR=' num2str(PR(j)), '  NSIVCB=', num2str(NSIVCB(j))]; 

TEXT3=['  AnMI=', num2str(AnMI(j)), '  MI=' num2str(MI(j)),  '  OldMI=', num2str(OldMI(j)),...
    '  RAD=', num2str(RAD(j)), '  LAD=', num2str(LAD(j)), '  LPR=', num2str(LPR(j)),...
    '  LQT=', num2str(LQT(j)), '  SA=', num2str(SA(j))];

TEXT4=['LQRSV=' num2str(LQRSV(j)), '  TAb=', num2str(TAb(j)), '  TInv=' , num2str(TInv(j)),...
    '  SB=', num2str(SB(j)), '  QAb=', num2str(QAb(j))];
    
fprintf('Output_1: %s\n',TEXT1);
fprintf('Output_2: %s\n',TEXT2);
fprintf('Output_3: %s\n',TEXT3);
fprintf('Output_4: %s\n',TEXT4);

else VF(j)=1; end
end

CRBBB(j)=0;

SVPB(j)=0;
VEB(j)=0;
% HDIAGN_star ={ 'IAVB' ;   'AF'    ; 'AFL'     ; 'Barady'  ; 'CRBBB'    ; 'IRBBB'  ; 'LAnFB'   ; 'LAD'     ; 'LBBB' ; ...
%                'LQRSV' ; 'NSIVCB' ; 'PR'      ;  'PAC'    ;  'PVC'     ;  'LPR'    ;  'LQT'   ; 'QaB'     ;  'RAD' ; ...
%                'RBBB'  ;  'SA'    ; 'SB'      ; 'SNR'     ; 'STach'    ; 'SVPB'   ; 'TAb'     ; 'Tinv'    ; 'VEB' ; '*Other*' } ;

my_label=[          I_AVB(j) AF(j) AFL(j) Brady(j) CRBBB(j) IRBBB(j) LAnFB(j) LAD(j) LBBB(j) ];
my_label=[my_label  LQRSV(j) NSIVCB(j) PR(j) PAC(j) PVC(j)  LPR(j) LQT(j) QAb(j) RAD(j) ];
my_label=[ my_label  RBBB(j) SA(j) SB(j) SNR(j) STach(j) SVPB(j) TAb(j) TInv(j) VEB(j)];

%T=['IAVB  AF AFL Brady IRBBB LAnFB LAD LBBB LQRSV NSIVCB PR PAC PVC LQT LPR QAb RAD RBBB SA SB SNR STach SVPB TAb TInv VEB']
T8=['       AF      Brady     LAnFB     LBBB     NSIVCB      PAC        LPR      QAb      RBBB       SB       STach     TAb       VEB'];
T9=['IAVB       AFL      IRBBB      LAD      LQRSV      PR        PVC       LQT       RAD       SA        SNR      SVPB      TInv    '];

fprintf('%s\n%s\n',T8,T9);
fprintf('%3.0f  ',my_label);fprintf('\n');

 %===================================   END program of classification ==================
% DIAGN_star=[270492004 	164889003 	164890007 	426627000 	713427006 	713426002 	445118002   39732003 	164909002 ...
%             251146004 	698252002 	10370003 	284470004 	427172004  164947007 	111975006 	164917005 	47665007  ...
%             59118001 	427393009 	426177001   426783006 	427084000 	63593006 	164934002 	59931005 	17338001 9999999 ];
% for i=1:numel(classes)
%     classes_2(i)=str2num(classes{i});
% end
        


function plot(a,b,c,d,e,f,g,h,i,j,k,l,m,k1,k2,k3,k4,k5,k6,k7,k8,k9,k10)
end
function subplot(a,b,c,d)
end
function figure(a,b,c)
end
function title(a,b,c,d)
end
function clf
end
function hold (a,b,c)
end
function axis(a,b,c,d,e,f)
end
function view(a,b,c,d,e,f,g,h)
end
function xlabel(a,b,c,d,f,g)
end
function ylabel(a,b,c,d,e,f,g)
end

