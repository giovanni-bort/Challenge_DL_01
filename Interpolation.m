function [Y]=Interpolation(X,Hz_init,Hz_final)

%Hz_int=1500
X_init=1:length(X);
X_final=Hz_init/Hz_final:Hz_init/Hz_final:length(X);
Y = interp1(X_init,X,X_final,'spline');
end
