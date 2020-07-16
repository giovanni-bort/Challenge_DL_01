%	---	F U N C T I O N		----
%	Searching a flat segment to the left
%	The amplitude criteria (uV - fixed) is addaptive
%	If a flat is not found in the interval -> Crit=Crit+1 uV and start again

% Entrance:
% Lead - a lead to work on
% From - from this point
% To - to this point
% Crit - amplitude criteria
% Flat - No of points where Up or Down criteria must be fullfill

% Exit:
% Pnt - beg point where flat segment found
%     = 0 - no flat segment	

function[Pnt,Crit1]=Flat_l(Lead,Crit,From,To,Flat);
To=max(Flat+1+1,To);    %******** MODIFIED 16.04.2020   Gio
Crit1=abs(Crit);
Pnt=0; a=0; 
while Pnt==0;
    if a<2; a=a+1;			%increase a
    else  Crit1=Crit1+0.003; %increase Crit with 3 uV
        if Crit1>0.500; Pnt=floor((From-To)/2); end % Flat_l not found
    end
    
    for i=From:-1:To;
         d(1:Flat)=Lead(i:-1:i-Flat+1)-Lead(i-1:-1:i-Flat);
         if max(abs(d)) < Crit1 &...
                 abs(Lead(i)-Lead(i-Flat-1)) < a*Crit1 &...
                 abs(Lead(i)-Lead(i-floor(Flat/2))) < 1.4*a*Crit1 &...	%1.4
                 abs(Lead(i-Flat)-Lead(i-floor(Flat/2))) < 1.4*a*Crit1;	%1.4
             Pnt=i; break
         end
    end
end
             
