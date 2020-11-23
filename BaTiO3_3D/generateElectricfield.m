function [E_step,E_coordinate]=generateElectricfield(magnitude,step)
    E_temp1=(0:(magnitude/step):magnitude);
    E_temp2=(magnitude:(-magnitude/step):-magnitude);
    E_temp2(1)=[];
    E_temp3=(-magnitude:(magnitude/step):0);
    E_temp3(1)=[];
    E_step=[E_temp1,E_temp2,E_temp3,E_temp1(2:length(E_temp1))];
    E_step=E_step(2:length(E_step));
    E_coordinate=(-magnitude:(magnitude/step):magnitude);
end
    
    
    