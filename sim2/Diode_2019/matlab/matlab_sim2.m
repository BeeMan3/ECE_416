
% Starter Code for Simulation Assignment #2
% Created by Michael Hayes 

% Note:  Some values below will need to be filled in.


clear all;
close all;
%Constants
Na = [1E19, 5E17, 1E16, 5E15];                 %p-side doping (in cm^-3)
Nd = [1E15, 1E16, 5E17, 1E18]; 				   %n-side doping (in cm^-3)
A = 3*10^(-8);                         %Area (in cm2)   
T = 300;                                       %Temperature (K)
Tau = 15*10^-6;                                %Minority carrier lifetime (sec)
Taud= 100*10^-6; 
M = 0.5;                   
ni = 2.1*(10^6);                                    %Intrinsic carrier conc (cm^-3)
Kgaas = 12.9;
e0 = 8.85*10^-14;                              %Permittivity of free space (F/cm)
q = 1.6*10^-19;                                % C
V_bias = -5;                                   %Bias voltage
Ksi = 11.8;
k = 8.617*10^(-5);

un1=500+(8900/(1+(Nd(1)/(6*(10^(16)))))^(.394));
up1=20+(471.5/(1+(Na(1)/(1.48*(10^(17)))))^(.38));
un2=500+(8900/(1+(Nd(2)/(6*(10^(16)))))^(.394));
up2=20+(471.5/(1+(Na(2)/(1.48*(10^(17)))))^(.38));
un3=500+(8900/(1+(Nd(3)/(6*(10^(16)))))^(.394));
up3=20+(471.5/(1+(Na(3)/(1.48*(10^(17)))))^(.38));
un4=500+(8900/(1+(Nd(4)/(6*(10^(16)))))^(.394));
up4=20+(471.5/(1+(Na(4)/(1.48*(10^(17)))))^(.38));
%PART A%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%i)
Dn1=(.0259)*un1;
Dp1=(.0259)*up1;
Ln1=sqrt(Dn1*Tau);
Lp1=sqrt(Dp1*Tau);
I01=q*A*(((Dn1/Ln1)*((ni^2)/Na(1)))+((Dp1/Lp1)*((ni^2)/Nd(1))));

Vbi1=(.0259)*(log((Na(1)*Nd(1))/(ni^2)));

W1n=sqrt(((2*Ksi*e0)/q)*((Na(1)+Nd(1))/(Na(1)*Nd(1)))*(Vbi1-0));
W1v=sqrt(((2*Ksi*e0)/q)*((Na(1)+Nd(1))/(Na(1)*Nd(1)))*(Vbi1-V_bias));

Cj01=(Ksi*e0*A)/W1n;
Cj1=(Ksi*e0*A)/W1v;

%ii)
Dn2=(.0259)*un2;
Dp2=(.0259)*up2;
Ln2=sqrt(Dn2*Tau);
Lp2=sqrt(Dp2*Tau);
I02=q*A*(((Dn2/Ln2)*((ni^2)/Na(2)))+((Dp2/Lp2)*((ni^2)/Nd(2))));

Vbi2=(.0259)*(log((Na(2)*Nd(2))/(ni^2)));

W2n=sqrt(((2*Ksi*e0)/q)*((Na(2)+Nd(2))/(Na(2)*Nd(2)))*(Vbi2-0));
W2v=sqrt(((2*Ksi*e0)/q)*((Na(2)+Nd(2))/(Na(2)*Nd(2)))*(Vbi2-V_bias));

Cj02=(Ksi*e0*A)/W2n;
Cj2=(Ksi*e0*A)/W2v;

%iii
Dn3=(.0259)*un3;
Dp3=(.0259)*up3;
Ln3=sqrt(Dn3*Tau);
Lp3=sqrt(Dp3*Tau);
I03=q*A*(((Dn3/Ln3)*((ni^2)/Na(3)))+((Dp3/Lp3)*((ni^2)/Nd(3))));

Vbi3=(.0259)*(log((Na(3)*Nd(3))/(ni^2)));

W3n=sqrt(((2*Ksi*e0)/q)*((Na(3)+Nd(3))/(Na(3)*Nd(3)))*(Vbi3-0));
W3v=sqrt(((2*Ksi*e0)/q)*((Na(3)+Nd(3))/(Na(3)*Nd(3)))*(Vbi3-V_bias));

Cj03=(Ksi*e0*A)/W3n;
Cj3=(Ksi*e0*A)/W3v;

%iv
Dn4=(.0259)*un4;
Dp4=(.0259)*up4;
Ln4=sqrt(Dn4*Tau);
Lp4=sqrt(Dp4*Tau);
I04=q*A*(((Dn4/Ln4)*((ni^2)/Na(4)))+((Dp4/Lp4)*((ni^2)/Nd(4))));

Vbi4=(.0259)*(log((Na(4)*Nd(4))/(ni^2)));

W4n=sqrt(((2*Ksi*e0)/q)*((Na(4)+Nd(4))/(Na(4)*Nd(4)))*(Vbi4-0));
W4v=sqrt(((2*Ksi*e0)/q)*((Na(4)+Nd(4))/(Na(4)*Nd(4)))*(Vbi4-V_bias));

Cj04=(Ksi*e0*A)/W4n;
Cj4=(Ksi*e0*A)/W4v;


%PART B%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V_graph = linspace(-1.5,1.5);
I1 = zeros(length(I01), length(V_graph));
I2 = zeros(length(I02), length(V_graph));
I3 = zeros(length(I03), length(V_graph));
I4 = zeros(length(I04), length(V_graph));
for i = 1:length(I01)
    I1(i,:) = I01(i).*(exp((V_graph)./(0.0259))-1);
end
for j = 1:length(I02)
    I2(j,:) = I02(j).*(exp((V_graph)./(0.0259))-1);
end
for m = 1:length(I03)
    I3(m,:) = I03(m).*(exp((V_graph)./(0.0259))-1);
end
for n = 1:length(I04)
    I4(n,:) = I04(n).*(exp((V_graph)./(0.0259))-1);
end
figure(1)
semilogy(V_graph,abs(I1),V_graph,abs(I2),V_graph,abs(I3),V_graph,abs(I4)); grid;
xlabel('Voltage'); ylabel('Current');
legend('Na=1e19 Nd=1e15','Na=5e17 Nd=1e16','Na=1e16 Nd=5e17','Na=5e15 Nd=1e18','location','northwest')

%PART D%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Dn3=(.0259)*un3;
Dp3=(.0259)*up3;
Ln3=sqrt(Dn3*Taud);
Lp3=sqrt(Dp3*Taud);
I03d=q*A*(((Dn3/Ln3)*((ni^2)/Na(3)))+((Dp3/Lp3)*((ni^2)/Nd(3))));

I3 = zeros(length(I03), length(V_graph));
I3d = zeros(length(I03d), length(V_graph));

for m = 1:length(I03)
    I3(m,:) = I03(m).*(exp((V_graph)./(0.0259))-1);
end
for m = 1:length(I03d)
    I3d(m,:) = I03d(m).*(exp((V_graph)./(0.0259))-1);
end
figure(2)
semilogy(V_graph,abs(I3),V_graph,abs(I3d)); grid;
xlabel('Voltage'); ylabel('Current');
legend('tau=15e-6','tau=100e-6')