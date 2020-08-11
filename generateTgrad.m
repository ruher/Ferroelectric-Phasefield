function T=generateTgrad(originT,finalT,interval,Tcycle)
    Tlength=finalT-originT;
    Tinterval=Tlength/interval;
    if Tcycle==true
        T=[originT:Tinterval:finalT,finalT-Tinterval:-Tinterval:originT];
    else
        T=(originT:Tinterval:finalT);
    end
end