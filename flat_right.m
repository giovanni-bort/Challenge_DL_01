% Entrance:
% Lead - a lead to work on
% From - from this point
% To - to this point
% Crit - amplitude criteria
% Flat - No of points where Up or Down criteria must be fullfill

% Exit:
% Pnt - beg point where flat segment found
% Crit - amplitude criteria has been changed

function[Pnt,Crit]=flat_right(Lead,Crit,From,To,Flat);
To=min(numel(Lead)-Flat-1,To);                %********** MODIFIED 16.04.2020 GB
Pnt=0; a=0;
    while Pnt==0;
        if a<2; a=a+1; 			%increase a
        else Crit=Crit+0.002;	%increase Crit with 1 uV
        end
        
        for ii=From:To;
            d(1:Flat)=Lead(ii:ii+Flat-1)-Lead(ii+1:ii+Flat);
            if max(abs(d)) < Crit &...
                    abs(Lead(ii)-Lead(ii+Flat+1)) < 4*a*Crit &...
                    abs(Lead(ii)-Lead(ii+floor(Flat/2))) < 3*a*Crit &...
                    abs(Lead(ii+Flat)-Lead(ii+floor(Flat/2))) < 3*a*Crit &...
                    abs(Lead(ii))<0.2+Crit;
                Pnt=ii; break
            end 
        end
    end