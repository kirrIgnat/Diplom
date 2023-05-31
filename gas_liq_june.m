clear all;
Nx=500; Nx1=Nx+1; Nx2=Nx/2;
xbeg=0.; xend=10.;
h=(xend-xbeg)/Nx;
  x=xbeg:h:xend; % grid is made
 roa=zeros(1,Nx); rob=zeros(1,Nx); p=zeros(1,Nx); pinf=zeros(1,Nx);
 u=zeros(1,Nx); E=zeros(1,Nx); T=zeros(1,Nx);
 F=0.; Q=0.;
gma = 1.4;  gmb = 2.8;
Sc=1.; Pr=0.7; Tbeg = 308.15;
gma = 1.4; pinfa = 0.;
Cpa = 1004.5; Ra = Cpa*(gma-1.)/gma; %Air
gmb = 2.8; pinfb=8.5e8;
Cpb = 4186.; Rb = Cpb*(gmb-1.)/gmb; %Water
alpha = 0.5; beta = 0.45;
time=0.; tfin=0.002;

%Initial conditions E - total energy
%Air
for i=1:Nx2
p(i)=1.e9;
T(i)=Tbeg; rob(i)=0.;
roa(i)=(p(i)+pinfa)/(Ra*T(i));
E(i)=(p(i)+gma*pinfa)/(gma-1);
pinf(i)=pinfa;
end
%Water
for i=(Nx2+1):Nx
p(i)=1.e5;
T(i)=Tbeg; roa(i)=0.;
rob(i)=(p(i)+pinfb)/(Rb*T(i));
E(i)=(p(i)+gmb*pinfb)/(gmb-1);
pinf(i)=pinfb;
end
ro=roa+rob;
% Сглаживание в средней точке
 % i=Nx/2;
 % roa(i)=(roa(i)+roa(i+1))*0.5;
 % rob(i)=(rob(i)+rob(i+1))*0.5;
 % p(i)=(p(i)+p(i+1))*0.5;
 % T(i)=(T(i)+T(i+1))*0.5;
 % E(i)=(E(i)+E(i+1))*0.5;
 % ro(i)=roa(i)+rob(i);
%mail loop    
 while time<tfin   
 % Calculation dt
 for i=1:Nx
  E1(i)=sqrt((gma*roa(i)*Ra+gmb*rob(i)*Rb)*T(i)/ro(i))+abs(u(i)); 
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
    Tm=(T(i)+T(i-1))*0.5;
    cva=Ra/(gma-1.); cpa=cva+Ra;
    cvb=Rb/(gmb-1.); cpb=cvb+Rb;
    cvm=(rma*cva+rmb*cvb)/rm;
    cpm=(rma*cpa+rmb*cpb)/rm;
    gam=cpm/cvm;
    gamCs=gma*rma*Ra+gmb*rmb*Rb;
    cs=sqrt(gamCs*Tm/rm);
    tau=alpha*h/cs; 
    visc=0.01*tau*pm*Sc;
    visc=0.;
    cond=visc/(Pr*(gam-1));
    dprux=(p(i)+ro(i)*u(i)*u(i)-p(i-1)-ro(i-1)*u(i-1)*u(i-1))/h;
    w=tau/rm*(dprux-rm*F);
      % upwind
    Froa=rma*(um-w); Frob=rmb*(um-w); Fro=Froa+Frob;
    
     % добавка в ур-ие неразрывности
    dprux=(p(i)-p(i-1))/h+rm*um*(u(i)-u(i-1))/h;
    w=tau/rm*(dprux-rm*F);
    Froa=rma*(um-w); Frob=rmb*(um-w);
           % конец уменьшения w
    Froa=Froa-tau*um*(u(i)*roa(i)-u(i-1)*roa(i-1))/h;
    Frob=Frob-tau*um*(u(i)*rob(i)-u(i-1)*rob(i-1))/h; 
      %  Froa=Froa-tau*rm*um*(u(i)-u(i-1))/h;
     %  Frob=Frob-tau*um*rm*(u(i)-u(i-1))/h; 
          % конец добавки в ур-ие неразрывности   
          
    dpx=(p(i)-p(i-1))/h; dux=(u(i)-u(i-1))/h;
    Wx=tau*(rm*um*dux+dpx-rm*F);
    Ptau=(um*dpx+gam*pm*dux-(gam-1)*Q)*tau;
    Pxx=4./3.*visc*dux+um*Wx+Ptau;
    Fu=pm+um*Fro-Pxx;
%    dex=(p(i)/ro(i)-p(i-1)/ro(i-1))/(h*(gam-1.)); % gam????
    dex=(E(i)/ro(i)-E(i-1)/ro(i-1))/h-um*(u(i)-u(i-1))/h; 
    dmr=(1./ro(i)-1./ro(i-1))/h;
    qst=tau*rm*um*(um*dex+pm*um*dmr-Q/rm);
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
% New variables     
dts=dt/h;
roa1=roa1*dts + roa; rob1=rob1*dts + rob;
ro1 = roa1 + rob1;
for i=2:Nx-1
  drux=(ro(i+1)*u(i+1)-ro(i-1)*u(i-1))/(2.*h);
    RSu=F*(ro(i)-tau*drux);
  dprux=(p(i+1)+ro(i+1)*u(i+1)*u(i+1)-p(i-1)-ro(i-1)*u(i-1)*u(i-1))/(2.*h);  
  w=tau*(dprux-ro(i)*F);
    RSE=Q+F*(ro(i)*u(i)-w);
  u1(i)=u1(i)*dts + ro(i)*u(i)+ dt*RSu;
  E(i)=E1(i)*dts + E(i) + dt*RSE;   
 end
 for i=2:Nx-1
   roa(i)=roa1(i); rob(i)=rob1(i); ro(i)=roa(i)+rob(i);
   u(i)=u1(i)/ro(i);
   cvr=roa(i)*Ra/(gma-1.)+rob(i)*Rb/(gmb-1.);
   pinf(i)=(roa(i)*pinfa+rob(i)*pinfb)/ro(i);
   T(i)=(E(i)-0.5*ro(i)*u(i)*u(i)-pinf(i))/cvr;
   p(i)=T(i)*(Ra*roa(i)+Rb*rob(i))-pinf(i);
 end
 end
 xx=xbeg+h/2:h:xend-h/2;
 plot(xx,p,'b','LineWidth',2); %pause(5);
 %plot(xx,ro,'b',xx,roa,'r',xx,rob,'g'); pause(5);
 % xlabel('x'); ylabel('Pressure'); 
  % plot(xx,ro,'b'); %pause(5);
% plot(xx,ro,'b',xx,roa,'r',xx,rob,'g'); pause(5);
  xlabel('x'); ylabel('Density'); 
 % r_res = vertcat(xx, ro, roa, rob).';
  fil_r=fopen('rho_12s.dat','w');
  con=roa./ro;
for i=1:Nx
    fprintf(fil_r,...
    '%10.6f %10.6f %10.6f %10.6f %10.6f %10.1f %10.6f %10.6f\n',...
    x(i),ro(i),roa(i),rob(i),u(i),p(i),T(i),con(i));
end
fclose(fil_r);
