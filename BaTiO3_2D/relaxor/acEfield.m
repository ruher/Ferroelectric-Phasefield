function E=acEfield(waveform,n,Emax,sampling)
t=(2*pi/sampling:(2*pi/sampling):2*n*pi);
if strcmp(waveform,'sin')==1
    E=Emax*sin(t);
elseif strcmp(waveform,'triangle')==1
    Etemp1=(Emax/sampling:(4*Emax/sampling):Emax);
    Etemp2=((Emax-Emax/sampling):(-4*Emax/sampling):-Emax);
    Etemp3=(-Emax:(4*Emax/sampling):0);
    Etemp=[Etemp1,Etemp2,Etemp3];
    E=[];
    for i=1:n
        E=[E,Etemp];
    end
elseif strcmp(waveform,'square')==1
    Etemp=[Emax*ones(1,ceil(sampling/2)),-Emax*ones(1,ceil(sampling/2))];
    E=[];
    for i=1:n
        E=[E,Etemp];
    end
end