clear,clc,
%i = 32.01;
%for i=1:.0001:2
 %   for j = 31.5:.0001:33.5
%initial conditions
R0km = 222;
v0=32;
o0=0;
y1=32;

% these values work
% y1                  32.013000
% anglee3            -0.053404
% alt3                  90.289371
%R0km = 222;
%v0=1.385430;


%Given Parameters
RE = 6378; %km
RM = 1738; %km
Rs =  66183; %km
D = 384400; %km
Uedu = 1; %du
Ue = 3.986*10^5; %km
wm = 2.137*10^-3;
Um = 4.9028*10^3;

TU = sqrt(RE^3/(Ue));




%Arrivial at SOI conditions
R0 = (R0km+RE)/RE;
Ddu = D / RE; %converts d to du
Rsdu = Rs/RE;
em0 = (v0^2/2)- (Uedu/R0);
Hm0 = R0*v0*cos(deg2rad(o0));
r1 = sqrt(Ddu^2+Rsdu^2-2*Ddu*Rsdu*cos(deg2rad(y1)));
v1 = sqrt(2*(em0+Uedu/r1));
o1 = acos(Hm0/(r1*v1));
phaseangle1 = asin((Rsdu/r1)*sin(deg2rad(y1)));
p0 = Hm0^2/Uedu; %du
a0 = -Uedu/(2*em0); %du
e0 =sqrt((a0-p0)/a0);

ta1 = acos(((p0/r1)-1)/e0);
E1 = 2*atan((sqrt((1-e0)/(1+e0))*tan(ta1/2)));
tof1 = sqrt(a0^3/Uedu)*(E1-e0*sin(E1)); %TU
tof1secs = tof1 * TU;
wmdt= wm*tof1;
phaseangledeparture = ta1-phaseangle1-wm*tof1;

%%%
DuTu = RE/TU;
v1kms = v1*DuTu;
r2 = Rs;
%Um = 4.9028*10^3; %km sec
%finds velocity of the moon
emm = -Ue/(2*D);
Vm = sqrt(2*(emm+(Ue/D)));



Ue = 3.986012*10^5;
v2 = sqrt(v1kms^2+Vm^2 -2*Vm*v1kms*cos(o1-phaseangle1));
angle2 = asin((Vm/v2)*cos(deg2rad(y1))-(v1kms/v2)*cos(deg2rad(y1)+phaseangle1-o1));
angle2deg = angle2*360/(2*pi); % the formula above might be wrong
em2 = (v2^2)/2-Um/r2;
%hmmoon = r2*v2*cos(90-o1*360/(2*pi));
hm2 = r2*v2*sin(angle2);
p2 = hm2^2/Um;
a2 = -Um/(2*em2); %du
e2 =sqrt((a2-p2)/a2);
rp3 = p2/(1+e2);
vp3 = sqrt(2*(em2+Um/rp3));
%disp(rp);
%disp(vp);
alt = rp3 - RM;
%disp(alt)
ta2 = acos(((p2/r2)-1)/e2) * -1; %multiply by -1 because going to perigee
%disp(rad2deg(ta2));
F2 = 2*atanh((sqrt((e2-1)/(1+e2))*tan(ta2/2)));
tof2 = sqrt((-a2)^3/Um)*(e2*sinh(F2)-F2);
TotalTOF = tof1secs/3600 + -1*tof2/3600;


if((alt > 89.9 && alt<90.1) && angle2<0)
fprintf('\nv0                  %f\n',v0);
fprintf('y1                  %f\n',y1);
fprintf('anglee3            %f\n',angle2);
fprintf('alt3                  %f\n',alt);
end
    %end
%end;
% 
% fprintf('variable              output\n');
% fprintf('v0                    %f\n',v0);
% fprintf('phi0                  %f\n',o0);
% fprintf('lamda1                %f\n',y1);
% fprintf('Energy0               %f\n',em0);
% ta0=0;%go back and fix/find this
% fprintf('true anomoly          %f\n',rad2deg(ta0));
% fprintf('true anomoly1         %f\n',rad2deg(ta1));
% fprintf('wmoondt               %f\n', wmdt);
% fprintf('phaseangle at depart  %f\n', rad2deg(phaseangledeparture));
% fprintf('r1                    %f\n',r1);
% fprintf('v1                    %f\n',r1);
% fprintf('o1                    %f\n',rad2deg(o1));
% fprintf('phaseangle arrival    %f\n',rad2deg(phaseangle1));
% fprintf('aCL                   %f\n',a0);
% fprintf('eCL                   %f\n',e0);
% fprintf('v2                    %f\n',v2);
% fprintf('o2                    %f\n',90-deg2rad(angle2)); %go back and look at this?
% fprintf('angle e2              %f\n',rad2deg(angle2));
% fprintf('energy2               %f\n',em2);
% fprintf('a2                    %f\n',a2);
% fprintf('e2                    %f\n',e2);
% fprintf('rp3                   %f\n',rp3);
% fprintf('alt3                  %f\n',alt);
% fprintf('vp3                   %f\n',vp3);
% 
% fprintf('\n anglee3            %f\n',angle2);
% fprintf('alt3                  %f\n',alt);