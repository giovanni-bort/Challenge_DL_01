%	---	F U N C T I O N		----
% Angle calculation in 360 degr (the 4th quadrants)

function[ang]=quadr(dy,dx);
ang=atan(abs(dy/dx))*180/pi;				% angle 

if dy>=0 & dx>=0;			ang=ang;
elseif dy>0 & dx<=0;		ang=180-ang;
elseif dy<0 & dx<0;		ang=ang+180;
elseif dy<0 & dx>=0;		ang=360-ang;
end
