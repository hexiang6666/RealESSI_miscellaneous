%%% main function for SPEC %%%%


%%%%% add by usr to provide input arguments %%%%%%%
clc; 

clear;

g=10.0;

Ag=load('acc.txt');

Ag = Ag(:,2);

% Ag=Ag;

dt=2.902139723300933838e-02;

zet=5.0;

endp= 5.0; 

[T,Spa,Spv,Sd]=SPEC(dt,Ag,zet,g,endp); 

frequency=1.0./T;

Spa=Spa.*g;

% plot(frequency, Spa);
% 
% xlabel('Frequency [HZ]');
% 
% ylabel('Pseudo acceleration [m/s^2]');

save 'Period.txt' T  -ascii;
save 'frequency.txt' frequency -ascii;
save 'Spa.txt' Spa -ascii;






