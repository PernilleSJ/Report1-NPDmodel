function I=lightfunc(P,param)

int=cumsum(P.*param.dz)*param.k;

I=param.I0*exp(-param.Kbg*param.z-int);

% Make I a column vector:
I = I';
end 