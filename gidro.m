function [result] = gidro(alpha,beta, Nx)
Nx1=Nx+1; Nx2=Nx/2;
xbeg=0.; xend=1.;
h=(xend-xbeg)/Nx;
x=xbeg:h:xend; % grid is made
roa=zeros(1,Nx); rob=zeros(1,Nx); p=zeros(1,Nx);
u=zeros(1,Nx); E=zeros(1,Nx); T=zeros(1,Nx);
F=0.; Q=0.;
gma = 5./3.;  gmb = 1.4; %gma = 1.4;
Sc=1.; Pr=0.7;
Runiv = 8314.;
Ra = Runiv/4.003; %Helium
Rb = Runiv/28.96; %Air

time=0.; tfin=0.0002;

%Initial conditions E - total energy
%Helium
for i=1:Nx2
    roa(i)=14.54903; rob(i)=0.;
    p(i)=1.943e7;
    T(i)=p(i)/(Ra*roa(i));
    E(i)=p(i)/(gma-1);
end
%Air
for i=(Nx2+1):Nx
    roa(i)=0.; rob(i)=1.16355;
    p(i)=100000.;
    T(i)=p(i)/(Rb*rob(i));
    E(i)=p(i)/(gmb-1);
end
ro=roa+rob;

%i=Nx/2;
%roa(i)=(roa(i)+roa(i+1))*0.5;
%rob(i)=(rob(i)+rob(i+1))*0.5;
%  p(i)=(p(i)+p(i+1))*0.5;
%  T(i)=(T(i)+T(i+1))*0.5;
%  E(i)=(E(i)+E(i+1))*0.5;
%  ro(i)=roa(i)+rob(i);

while time<tfin   
 % Calculation dt
    for i=1:Nx
        E1(i)=sqrt((gma*roa(i)+gmb*rob(i))*p(i)/(ro(i)*ro(i)))+abs(u(i)); 
    end
    dt=beta*h/max(E1);
    time=time+dt;
    for i=2:Nx
        roa1(i)=0.; rob1(i)=0.; u1(i)=0.; E1(i)=0.;
    end
% Fluxes 
    for i=2:Nx
       rm=(ro(i)+ro(i-1))*0.5; 
    % upwind
       rmaf=roa(i-1); 
       rma=(roa(i)+roa(i-1))*0.5;  
       rmbf=rob(i-1); 
       rmb=(rob(i)+rob(i-1))*0.5;
      % end upwind
        um=(u(i)+u(i-1))*0.5;
        pm=(p(i)+p(i-1))*0.5;
        Em=(E(i)+E(i-1))*0.5;
        cva=Ra/(gma-1.); cpa=cva+Ra;
        cvb=Rb/(gmb-1.); cpb=cvb+Rb;
        cvm=(rma*cva+rmb*cvb)/rm;
        cpm=(rma*cpa+rmb*cpb)/rm;
        gam=cpm/cvm;
        gamCs=(gma*rma*Ra+gmb*rmb*Rb)/(Ra*rma+Rb*rmb);
        cs=sqrt(gamCs*pm/rm);
    %  cs=sqrt(gam*pm/rm);
        tau = (alpha*h/cs); 
        visc=30*tau*pm*Sc;
       visc=0.;
        cond=visc/(Pr*(gam-1));
        dprux=(p(i)+ro(i)*u(i)*u(i)-p(i-1)-ro(i-1)*u(i-1)*u(i-1))/h;
        dpx=(p(i)-p(i-1))/h; dux=(u(i)-u(i-1))/h;
        wm=  (tau/rm*(um * rm * dux+dpx - rm*F));
        w=  wm+(tau/rm*(um * rm * dux+dpx - rm*F));

        %w=  (tau/rm*(um * rm * dux+dpx));
      % upwind
        Froa=rma*(um-w); Frob=rmbf*(um-w); Fro=Froa+Frob;
    
        Wmx=  (tau/rm*(um * rm * dux+dpx - rm*F));
        Wx= Wmx+(tau/rm*(um * rm * dux+dpx - rm*F));
        Ptau=(um*dpx+gam*pm*dux-(gam-1)*Q)*tau;
        Pxx=4./3.*visc*dux+um*Wx+Ptau;
        Fu=(pm+um*Fro-Pxx);
%    dex=(p(i)/ro(i)-p(i-1)/ro(i-1))/(h*(gam-1.)); % gam????
        dex=(E(i)/ro(i)-E(i-1)/ro(i-1))/h-um*(u(i)-u(i-1))/h; 
        dmr=(1./ro(i)-1./ro(i-1))/h;
        qst=(tau*rm*um*(um*dex+pm*um*dmr-Q/rm));
        FE=(Em+pm)*Fro/rm-cond*(T(i)-T(i-1))/h-qst-Pxx*um;
        roa1(i)=roa1(i)+Froa;
        roa1(i-1)=roa1(i-1)-Froa;
        rob1(i)=rob1(i)+Frob;
        rob1(i-1)=rob1(i-1)-Frob;
        u1(i)=u1(i)+Fu;
        u1(i-1)=u1(i-1)-Fu;
        E1(i)=E1(i)+FE;
        E1(i-1)=E1(i-1)-FE;
    end
    dts=dt/h;
    roa1=roa1*dts + roa; rob1=rob1*dts + rob;
    ro1 = roa1 + rob1;
    for i=2:Nx-1
        drux=(ro(i+1)*u(i+1)-ro(i-1)*u(i-1))/(2.*h);
        RSu=F*(ro(i)-tau*drux);
        dprux=(p(i+1)+ro(i+1)*u(i+1)*u(i+1)-p(i-1)-ro(i-1)*u(i-1)*u(i-1))/(2.*h);  
        wmm =(tau/rm*(um * rm * dux+dpx - rm*F));
        w=(tau/rm*(um * rm * dux+dpx - rm*F));
        RSE=Q+F*(ro(i)*u(i)-w);
        u1(i)=u1(i)*dts + ro(i)*u(i)+ dt*RSu;
        E(i)=E1(i)*dts + E(i) + dt*RSE;   
    end
    for i=2:Nx-1
        roa(i)=roa1(i); rob(i)=rob1(i); ro(i)=roa(i)+rob(i);
        u(i)=u1(i)/ro(i);
        cvr=roa(i)*Ra/(gma-1.)+rob(i)*Rb/(gmb-1.);
        T(i)=(E(i)-0.5*ro(i)*u(i)*u(i))/cvr;
        p(i)=T(i)*(Ra*roa(i)+Rb*rob(i));
    end
end
xx=xbeg+h/2:h:xend-h/2;
result = [xx; ro];
end

