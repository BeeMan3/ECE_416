%PART 1%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
k = 8.617e-5;

for ii=1:5;
    if ii == 1
        T = .000000000000000000000000001; %T is approxmately 0
    end
    if ii == 2
        T = 150;
    end
    if ii == 3
        T = 250;
    end
    if ii == 4
        T = 350;
    end
    if ii == 5
        T = 450;
    end
    
    kT = k*T;
    dE(ii,1)=-5*kT;
    for jj=1:101
        f(ii,jj)=1/(1+exp(dE(ii,jj)/kT)); %solve for fermi energy based on i and j 
        dE(ii,jj+1)=dE(ii,jj)+0.1*kT; %increment dE
    end
end
dE=dE(:,1:jj);

close
figure(1)
plot(dE',f');grid;
xlabel('E-EF(ev)'); ylabel('f(E)');
legend('T=0','T=150','T=250','T=350','T=450')
%PART 2%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NDref = 6.95e18;NAref=3.75e18;
unmin=88;upmin=54.3;
un0=1252;up0=407;

N = logspace(14,20); %set up logspace from e14 to e20
un = unmin + un0./(1+(N/NDref)); %solve for mobility of n type
up = upmin + up0./(1+(N/NAref)); %solve for mobility of p type


figure(2)
loglog(N,un,N,up);grid;
axis([1.0e14 1.0e20 1.0e1 1.0e4]);
xlabel('NA or ND (cm^3)');
ylabel('Mobility (cm^2/v*sec)');
text(1.0e15,1530,'Electrons');
text(1.0e15,530,'Holes');
text(1.0e18,2000,'SI,300K');

%PART 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ND1 = 10^(10);
ND2 = 10^(14)+10^(10);
ND3 = 10^(16)+10^(10);
ND4 = 10^(19)+10^(10);
NA1 = 10^(10);
NA2 = 10^(14)+10^(10);
NA3 = 10^(16)+10^(10);
NA4 = 10^(19)+10^(10);

q = 1.6e-19;
t = 0.01;
w = 0.01;
L = 0.0075;
unmin=88;upmin=54.3;
un0=1252;up0=407;

%solve for mobilities
un1 = unmin + un0./(1+(ND1/NDref));
up1 = upmin + up0./(1+(NA1/NAref));
un2 = unmin + un0./(1+(ND2/NDref));
up2 = upmin + up0./(1+(NA2/NAref));
un3 = unmin + un0./(1+(ND3/NDref));
up3 = upmin + up0./(1+(NA3/NAref));
un4 = unmin + un0./(1+(ND4/NDref));
up4 = upmin + up0./(1+(NA4/NDref));


ronarr1 = (q*un1*ND1);
ron1 = [];
ronarr2 = (q*un2*ND2);
ron2 = [];
ronarr3 = (q*un3*ND3);
ron3 = [];
ronarr4 = (q*un4*ND4);
ron4 = [];
ronarr5 = (q*up1*NA1);
ron5 = [];
ronarr6 = (q*up2*NA2);
ron6 = [];
ronarr7 = (q*up3*NA3);
ron7 = [];
ronarr8 = (q*up4*NA4);
ron8 = [];

%solve for ro 
for i = 1:length(ronarr1)
    ron1(i) = 1/ronarr1(i);
end
for i = 1:length(ronarr2)
    ron2(i) = 1/ronarr2(i);
end
for i = 1:length(ronarr3)
    ron3(i) = 1/ronarr3(i);
end
for i = 1:length(ronarr4)
    ron4(i) = 1/ronarr4(i);
end
for i = 1:length(ronarr5)
    ron5(i) = 1/ronarr5(i);
end
for i = 1:length(ronarr6)
    ron6(i) = 1/ronarr6(i);
end
for i = 1:length(ronarr7)
    ron7(i) = 1/ronarr7(i);
end
for i = 1:length(ronarr8)
    ron8(i) = 1/ronarr8(i);
end

V = 0:0.001:1;

%find resistance and current values
R1 = ((ron1/t)*(L/w));
I1 = V/R1;

R2 = ((ron2/t)*(L/w));
I2 = V/R2;

R3 = ((ron3/t)*(L/w));
I3 = V/R3

R4 = ((ron4/t)*(L/w));
I4 = V/R4;

R5 = ((ron5/t)*(L/w));
I5 = V/R5;

R6 = ((ron6/t)*(L/w));
I6 = V/R6;

R7 = ((ron7/t)*(L/w));
I7 = V/R7;

R8 = ((ron8/t)*(L/w));
I8 = V/R8;

figure(3)
loglog(V,I1,V,I2,V,I3,V,I4,V,I5,V,I6,V,I7,V,I8);grid; %plot all curves on loglog graph
axis([0 1.0 0 1e7]);
xlabel('Voltage')
ylabel('Current')
legend('Nd = 1e10','Nd = 1e14','Nd = 1e16','Nd = 1e19','Na = 1e10','Na = 1e14','Na = 1e16','Na = 1e19')


