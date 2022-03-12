function I=lightfuncseason(t,P,param)

int=cumsum(P.*param.dz)*param.k;

%I=param.I0*exp(-param.Kbg*param.z-int);
%I=param.I0*exp(-param.Kbg*param.z-int')*(1+cos(2*pi.*t./365));
I=(param.I0.*(1-cos((2*pi.*t)./365))).*exp(-param.Kbg.*param.z-int);

% Make I a column vector:
I = I';
end 