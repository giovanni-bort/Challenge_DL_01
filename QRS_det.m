% QRS Detection 
% Created by Ivaylo Christov
%
% Publications:
% 1. Christov II (2004) Real time electrocardiogram QRS detection using combined 
% adaptive threshold, Biomedical Engineering Online, 3, 28, 
% http://www.biomedical-engineering-online.com/content/3/1/28
%
% 2. Christov I (2007) Assessment of the performance of the adaptive
% thresholding algorithm for QRS detection with the use of AHA database.
% Bioautomation, 6, pp. 27-37, 
% http://www.clbme.bas.bg/bioautomation/2007/vol_6.1/files/6_4.1.pdf 



function[QRS]=QRS_det(d1,d2);
Hz=250;

d1=d1-d1(1);      d2=d2-d2(1);


% Moving averaging filter with a first zero at the mains frequency 50Hz (60 Hz)
% (Sum of all samples in 20 ms interval) divided by (No of samples)
ECG1(3:length(d1)-2)=(d1(1:length(d1)-4)+d1(2:length(d1)-3)+d1(3:length(d1)-2)+...
    d1(4:length(d1)-1)+d1(5:length(d1)-0))/5;
ECG2(3:length(d2)-2)=(d2(1:length(d2)-4)+d2(2:length(d2)-3)+d2(3:length(d2)-2)+...
    d2(4:length(d2)-1)+d2(5:length(d2)-0))/5;

d1=ECG1; d2=ECG2;

% Moving averaging filtration in 28 ms interval (7 samples)
ECG1(4:length(d1)-3)=(d1(1:length(d1)-6)+d1(2:length(d1)-5)+d1(3:length(d1)-4)+...
    d1(4:length(d1)-3)+d1(5:length(d1)-2)+d1(6:length(d1)-1)+d1(7:length(d1)-0))/7;
ECG2(4:length(d2)-3)=(d2(1:length(d2)-6)+d2(2:length(d2)-5)+d2(3:length(d2)-4)+...
    d2(4:length(d2)-3)+d2(5:length(d2)-2)+d2(6:length(d2)-1)+d2(7:length(d2)-0))/7;


% Baseline wander suppression filter
% Forward high-pass recursive filter useing the formula:
%                       Y(i)=C1*(X(i)-X(i-1)) + C2*Y(i-1)
T=1/Hz;        % [s] - sampling period
Fc=0.625;            % [Hz]
C1=1/[1+tan(Fc*pi*T)];  C2=[1-tan(Fc*pi*T)]/[1+tan(Fc*pi*T)];
b=[C1 -C1]; a=[1 -C2];
ECG1_show=filter(b,a,ECG1);     ECG2_show=filter(b,a,ECG2);
Fc=2.2;               % [Hz]
C1=1/[1+tan(Fc*pi*T)];  C2=[1-tan(Fc*pi*T)]/[1+tan(Fc*pi*T)];
b=[C1 -C1]; a=[1 -C2];
ECG1=filter(b,a,ECG1);     ECG2=filter(b,a,ECG2);
    
% -------------------------------------------------------------
%                       QRS triggering from 2 channel ECG
%           3 trasholds: M - Steep slope
%                        F - integrated covering of the abs(Derivative)
%                        R - expected beat trap
%--------------------------------------------------------------
    
Y(3:length(ECG1)-2) = abs(ECG1(5:length(ECG1))-ECG1(1:length(ECG1)-4))...
    + abs(ECG2(5:length(ECG2))-ECG2(1:length(ECG2)-4));
Y(1)=Y(3); Y(2)=Y(3); Y(length(ECG1)-1)=Y(length(Y)); Y(length(ECG1))=Y(length(Y));
    
% Moving averaging filtration in 40 ms interval (10 samples for 250Hz sampling rate)
% The No of samples should be odd, i.e. 11 samples
 D(6:length(Y)-5)=(Y(1:length(Y)-10)+Y(2:length(Y)-9)+Y(3:length(Y)-8)+...
        Y(4:length(Y)-7)+Y(5:length(Y)-6)+Y(6:length(Y)-5)+Y(7:length(Y)-4)+...
        Y(8:length(Y)-3)+Y(9:length(Y)-2)+Y(10:length(Y)-1)+Y(11:length(Y)-0))/11;
    D(length(D)+1:length(D)+5)=0;

%               Initialization
[Max,ma]=max(D(Hz:3*Hz)); M=0.6*Max;  % ititial value of M between 1s-3s
if length(find(D(Hz+ma+.4*Hz:Hz+ma+.4*Hz+1.6*Hz)>M))==0
    M=0.5*Max;                          % we have selected ES with a big slew rate
end                                     % and should correct it
A=find(D(Hz:5*Hz)>M);
if isempty(A)==1; QRS=Hz;
else QRS=A(1)+Hz;
end
Mnew=M; Mold=M;


%       Buffers Initialization
RR=[.7*Hz .7*Hz .7*Hz .7*Hz .7*Hz];        % 800ms for the RR intervals
MM=[M M M M M];             % last 5 values
pos1=1:5; pos2=1:5; neg1=1:5;  neg2=1:5;

M(1:length(D))=M;           % to speed-up
R(1:length(D))=0;           % to speed-up
QRS(1:1000)=QRS;            % to speed-up

F=mean(D(QRS(3)-round(Hz*.35):QRS(3)));  % 350 ms
F(1:length(D))=F;           % to speed-up
MFR=M+F;
MaxD=0; pMaxD=QRS(3);       % Initialization of max Deriv & pointer for max Deriv


%       Counters Initialization
j=4;                        % QRS No


%               Algorithm's body
for k=QRS(3):length(D)-round(.300*Hz);                      % stop 300 ms before the length
    F(k-1)=F(k-2) + (max(D(k-round(Hz*.05):k-1))...                % integrated cover in 50ms
        - max(D(k-round(Hz*.35):k-1-round(Hz*.35)+round(Hz*.05))))/round(Hz*.6);
    MFR(k-1)=F(k-1)*0.4+M(k-1)+R(k-1);
    if MFR(k-1)<0.25*Mold;  MFR(k-1)=0.25*Mold; end  % do not allow MFR -> 0
    
      
    if k-QRS(j-1)<round(.250*Hz);
        M(k)=M(k-1);           % 250 (now 200ms) after gurrent QRS - no decrease of M, no test
    else
        if D(k)>MFR(k-1);
             
            %                                  QRS !
            % Correct the location of QRS with respect to the expected place
            % & Existence of bigger Max of D in 250 ms interval
        
            [Max,ma]=max(D(k:k+round(.250*Hz)));        % max in 250 ms interval
            Mnew=0.70*Max;
            if Mnew>2*mean(MM); Mnew=1.1*Mold; end;
            if k-QRS(j-1)-mean(RR)<0; % expected place ?
                QRS(j)=k+ma;
            else QRS(j)=k;
            end
                      
            RR(j-5*floor(j/5)+1)=QRS(j-1)-QRS(j-2);
            MM(j-5*floor(j/5)+1)=Mnew;
            MaxD=0;    % Initialization of max Deriv
            Mold=mean(MM); M(k)=mean(MM);
            j=j+1;
            
            % Back Correction of Extrasystola which is not detected:
            % A3001, A3010, A4002, A4004, A4007, A4009, A5001, A5007, A6005, A7001
            t1=QRS(j-2)-QRS(j-3);  t2=QRS(j-1)-QRS(j-2); 
            if (t1> mean(RR) | mean(RR)-t1<0.12*mean(RR)) & abs(t2-2*mean(RR))<0.5*mean(RR);
                
                % Search for sharpest peak in last RR buffer
                [Po1,P1]=max((ECG1(QRS(j-2)+round(0.2*Hz)+2:QRS(j-1)-round(0.2*Hz)-2) -...
                    ECG1(QRS(j-2)+round(0.2*Hz):QRS(j-1)-round(0.2*Hz)-4)).*...
                    (ECG1(QRS(j-2)+round(0.2*Hz)+2:QRS(j-1)-round(0.2*Hz)-2) -...
                    ECG1(QRS(j-2)+round(0.2*Hz)+4:QRS(j-1)-round(0.2*Hz))));
                P1=P1+QRS(j-2)+round(0.2*Hz)+2-1;
                    
                [Po2,P2]=max((ECG2(QRS(j-2)+round(0.2*Hz)+2:QRS(j-1)-round(0.2*Hz)-2) -...
                    ECG2(QRS(j-2)+round(0.2*Hz):QRS(j-1)-round(0.2*Hz)-4)).*...
                    (ECG2(QRS(j-2)+round(0.2*Hz)+2:QRS(j-1)-round(0.2*Hz)-2) -...
                    ECG2(QRS(j-2)+round(0.2*Hz)+4:QRS(j-1)-round(0.2*Hz))));
                P2=P2+QRS(j-2)+round(0.2*Hz)+2-1;
                if Po1>Po2; Po=Po1; P=P1; else Po=Po2; P=P2; end
                
                % Search for max in last RR buffer
                [MaxD pMaxD]=max(D(QRS(j-2)+round(0.2*Hz):QRS(j-1)-round(0.2*Hz)));
                pMaxD=pMaxD+QRS(j-2)+round(0.2*Hz);
                               
                if Po>0.004 & 3*MaxD>mean(MM); QRS(j)=QRS(j-1); QRS(j-1)=P; j=j+1; end
            end
      
        else                        % Outside QRS
            % Old M(k)=M(k-1)-Mnew*0.6/((1.2-0.2)*Hz)
            % Now M(k)=M(k-1)-Mnew*5/((1.2-0.2)*Hz) for HR >100
            if k-QRS(j-1)<1.000*Hz; M(k)=M(k-1)-Mnew*0.6/((1.2-0.2)*Hz);
            else M(k)=M(k-1);
            end
    
            % Old R(k)=R(k-1)-Mnew*.4/((1.2-0.2)*Hz)
            % Now R(k)=R(k-1)-Mnew*5/((1.2-0.2)*Hz) for HR >100
            if k-QRS(j-1)>mean(RR)/3 & k-QRS(j-1)<mean(RR) ; R(k)=R(k-1)-Mnew*0.4/((1.2-0.2)*Hz);
            else R(k)=R(k-1);
            end
    
            if MaxD<D(k); pMaxD=k; end   % Max of Deriv
            
        end   % if D(k)>MFR(k-1);
    end   % if k-QRS(j-1)<round(.200*Hz);
end   % for k=QRS(3):length(D)-round(.300*Hz);
           

% Find 1 or 2 QRSs at the begining
A=find(D(1:QRS(3)-.100*Hz)>M(1));
if isempty(A)==0    %not empty
    QRS(2)=A(1);
    for i=1:length(A)-1;
        if A(i+1)-A(i)>1; QRS(1)=QRS(2); QRS(2)=A(i+1);
        end
    end
end

if QRS(2)>QRS(1); QRS=QRS(1:j-1);
else QRS=QRS(2:j-1);
end
QRS1=QRS/Hz;                  % discretes->s


