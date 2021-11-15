%*************************************************************************%
% This script was developed by Wenyu Jiang and Cheng Lin                  %
% Department of Civil Engineering                                         %
% University of Victoria, Canada                                          %
% Email: wenyujiang@uvic.ca; chenglin918@uvic.ca                          %
% Date: 2021-06-26                                                        %
%*************************************************************************%
%*************************************************************************%
% The code is for analyzing a leading row pile in a pile group under local%
% scour and lateral loaded conditions.                                    %
% Soil-pile interactions are simulated using Matlock's p-y curve.         %                                                                
% Note: The script is for the following boundary conditions:              %
% 1. given pile-top lateral displacement;                                 %
% 2. given pile-top bending moment;                                       %
% 3. pile-tip bending moment is zero (flexible pile);                     %
% 4. pile-tip shear force is zero (flexible pile).                        %
%*************************************************************************%
close all
clear all
clc
tic
%pile parameters
Ncol=5;%no. of piles on x-dir.
Nrow=3;%no. of piles on y-dir. 
D=0.324;%pile diameter(m)
ax=3.92*D;%pile spacing on x-dir.(m)
ay=3.29*D;%pile spacing on y-dir.(m)
R=28600.0;%flexural rigidity of pile,R=EI(kN-m2)
L=12.07;%pile length(m)
Lco=0.495;%unsupported length(m)
Leo=L-Lco;%embedded length(m)
Q=0;%axial load at pile head(kN)
fm=0.6;%p-multiplier to account for geotechnical group effect
SSIeta=1.0;%correction factor to account for soil-pile interaction
%scour-hole parameters
Sd=6.0593*D;%scour-hole depth(m)
Swb=0.0*D;%scour-hole bottom width(m)
Sq=26.6;%scour-hole slope angle(deg)
Swt=Swb+Sd*cotd(Sq);%scour-hole top width(m)
Le=Leo-Sd;%post-scour embedded length(m)
Lc=L-Le;%post-scour unsupported length(m)
%soil parameters:
cu=43.565;%undrained shear strength(kPa)
Jclay=0.5;%coefficient in Matlock p-y curve (typically taken as 0.5)
uw=9.05;%effective unit weight(kN/m3)
e50=0.01;%strain at 50% of stress(empirical value, typically taken as 0.01)
y50=2.5*e50*D;%deflection at 50% of ultimate resistance(m)
kint=50000;%a guess value of the initital slope of a linear p-y relation for elastic solution(kN/m2=kPa)
%loading parameters: 
%set two expected deflection values(1 inch and 89 mm):
yt_1=0.0254;%1st deflection applied at pile head(m)
yt_2=0.089;%2nd deflection applied at pile head(m)
nstep_1=7;%number of loading steps to reach yt_1
nstep_2=13;%number of loading steps to reach yt_2
nstep=nstep_1+nstep_2;%total number of loading steps
dltyt_1=yt_1/(nstep_1-1);%incremental deflection to reach yt_1 from 0(m)
dltyt_2=(yt_2-yt_1)/(nstep_2-1);%incremental deflection to reach yt_2 from yt_1(m)
yt_vector=[linspace(dltyt_1,yt_1,nstep_1),linspace(yt_1+dltyt_2,yt_2,nstep_2)]';%(m)
%accrodingly, set two expected bending moment values:
Mt_1=0.0;%1st bending moment applied at pile head(kN-m)
Mt_2=0.0;%2nd bending moment applied at pile head(kN-m)
dltMt_1=Mt_1/(nstep_1-1);%incremental moment to reach Mt_1 from 0(kN-m)
dltMt_2=(Mt_2-Mt_1)/(nstep_2-1);%incremental moment to reach Mt_2 from Mt_1(kN-m)
Mt_vector=[linspace(dltMt_1,Mt_1,nstep_1),linspace(Mt_1+dltMt_2,Mt_2,nstep_2)]';%(kN-m)
%mesh parameters:
esdsr=0.05;%desired element size(m)
n=ceil(L/esdsr);%number of elements
h=L/n;%element size(m)
zpt=linspace(0,L,n+1)';%depth below pile top(m)
rawzblwsc=zpt-Lc;%depth below scour bottom(m)
zblwsc=rawzblwsc(rawzblwsc>=0);%(m)
Nndblwsc=length(zblwsc);%number of embedded nodes below scour bottom
for page=1:Nndblwsc
    if zblwsc(page)==0.0
        zblwsc(page)=0.0001;
    else
    end
end
ze_wg=zeros(Nndblwsc,1);%equivalent depth from wedge model(m)
for page=1:Nndblwsc
    [ze_wg(page)]=zeWedge(zblwsc(page),Sq,Sd,Swb,D,uw,cu);%(m)
end
ie1=(n+3)-Nndblwsc+1;%node tag of the 1st embedded node in a n+5 vector
zeTail_wg=[linspace(0,0,ie1-1)';ze_wg;0;0];%tailored vector(length=n+5)
yeNOtailor_allstep=zeros(n+5,nstep);%elastic lateral displacement including four imaginary nodes' results(m)
ye_allstep=zeros(n+1,nstep);%elastic lateral displacement(m)
Me_allstep=zeros(n+1,nstep);%elastic bending moment(kN-m)
Ve_allstep=zeros(n+1,nstep);%elastic shear force(kN)
Se_allstep=zeros(n+1,nstep);%elastic rotation(deg)
pe_allstep=zeros(n+1,nstep);%elastic soil resistance(kN/m)
for loadstep=1:nstep
    [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=SCeClayfmYtMt(zeTail_wg,fm,n,R,Q,h,Lco,kint,yt_vector(loadstep),Mt_vector(loadstep));
    yeNOtailor_allstep(:,loadstep)=outputye_notailor;%(m)
    ye_allstep(:,loadstep)=outputye;%(m)
    Me_allstep(:,loadstep)=outputMe;%(kN-m)
    Ve_allstep(:,loadstep)=outputVe;%(kN)
    Se_allstep(:,loadstep)=outputSe;%(deg)
    pe_allstep(:,loadstep)=outputpe;%(kN/m)
end
xnorm=outputx./D;%normalized depth below pre-scour ground(xD)
y_allstep_wg=zeros(n+1,nstep);%lateral displacement(m)
M_allstep_wg=zeros(n+1,nstep);%bending moment(kN-m)
V_allstep_wg=zeros(n+1,nstep);%shear force(kN)
S_allstep_wg=zeros(n+1,nstep);%rotation(deg)
p_allstep_wg=zeros(n+1,nstep);%soil resistance(kN/m)
k_allstep_wg=zeros(n+1,nstep);%k=p/y (kPa)
pu_allstep_wg=zeros(n+1,1);%ultimate soil resistance(kN/m)
for loadstep=1:nstep
    [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=...
        SCMatlockfmYtMt(zeTail_wg,fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt_vector(loadstep),Mt_vector(loadstep),yeNOtailor_allstep(:,loadstep));
    y_allstep_wg(:,loadstep)=outputy;%(m)
    M_allstep_wg(:,loadstep)=outputM;%(kN-m)
    V_allstep_wg(:,loadstep)=outputV;%(kN)
    S_allstep_wg(:,loadstep)=outputS;%(deg)
    p_allstep_wg(:,loadstep)=outputp;%(kN/m)
    pu_allstep_wg(:,1)=outputpu;%(kN/m)
    k_allstep_wg(:,loadstep)=outputk;%(kPa)
end
PTOPdisp_wg=[0,y_allstep_wg(1,:)]';%pile-top deflection(m)
PTOPshear_wg=[0,V_allstep_wg(1,:)]';%pile-top shear force(kN)
toc
%plot the computed results:
figure(1)
subplot(2,4,1)
box on
grid on
plot(y_allstep_wg(:,end),xnorm,'-r')
axis ij
xlabel('y(m)')
ylabel('Depth below pre-scour ground(xD)')
title('Deflection')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,2)
box on
grid on
plot(S_allstep_wg(:,end),xnorm,'-r')
axis ij
xlabel('\theta(^{\circ})')
ylabel('Depth below pre-scour ground(xD)')
title('Rotation')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,3)
box on
grid on
plot(M_allstep_wg(:,end),xnorm,'-r')
axis ij
xlabel('M(kN-m)')
ylabel('Depth below pre-scour ground(xD)')
title('Bending Moment')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,4)
box on
grid on
plot(V_allstep_wg(:,end),xnorm,'-r')
hold off
axis ij
xlabel('V(kN)')
ylabel('Depth below pre-scour ground(xD)')
title('Shear force')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,5)
box on
grid on
plot(p_allstep_wg(:,end),xnorm,'-r')
hold on
plot(pu_allstep_wg(:,1),xnorm,'-b')
hold off
axis ij
legend('p','p_{u}')
xlabel('p(kN/m)')
ylabel('Depth below pre-scour ground(xD)')
title('Soil Resistance')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,6)
box on
grid on
plot(p_allstep_wg(:,end)./y_allstep_wg(:,end),xnorm,'-r')
axis ij
xlabel('k=p/y(kPa)')
ylabel('Depth below pre-scour ground(xD)')
title('Stiffness k=p/y')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,7)
box on
grid on
plot(PTOPdisp_wg,PTOPshear_wg,'-o')
xlabel('y(m)')
ylabel('V(kN)')
title('Lateral load-deflection curve at pile head')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,8)
box on
grid on
plot([0,abs(M_allstep_wg(1,:))./R]',[0,abs(M_allstep_wg(1,:))]','-o')
xlabel('Curvature')
ylabel('Bending moment(kN-m)')
title('Bending moment-curvature curve at pile head')
set(gca,'FontName','Times','FontSize',11)
%===function #1 ===
function [ze]=zeWedge(z,Sq,Sd,Swb,D,uw,cu)
%z=depth below scour-hole bottom(m)
%Sq=scour-hole slope angle(deg)
%Sd=scour-hole depth(m)
%Swb=distance from periphery of scour-hole bottom to the periphery of the
%piles in the leading row of a pile group(m) 
%D=pile diameter(m)
%uw=submerged unit weigth of soil(kN/m3)
%cu=undrained shear strength of clay(kPa)
TTA=tand(Sq)/(1.0-tand(Sq));
F0=(0.5*uw*D+sqrt(2.0)*cu)*(z^2)+2.0*cu*D*z;
F1=(0.5*uw*D+sqrt(2.0)*cu)*(z^2+TTA*((z-Swb)^2))+2.0*cu*D*(z+TTA*(z-Swb)); 
F2=(0.5*uw*D+sqrt(2.0)*cu)*((z+Sd)^2-Sd*(2.0*Swb+Sd*cotd(Sq)))+2.0*cu*D*(z+Sd);
if Sq<45.0 
    if z>=0 && z<=Swb 
        Fsc=F0;
    elseif z>Swb && z<=Swb+Sd/TTA 
        Fsc=F1;
    elseif z>Swb+Sd/TTA 
        Fsc=F2;
    end
else 
    if z>=0 && z<=Swb 
        Fsc=F0;
    elseif z>Swb 
        Fsc=F2;
    end
end
a_ze=0.5*uw*D+sqrt(2.0)*cu;
b_ze=2.0*cu*D;
c_ze=-Fsc;
Delta_ze=b_ze^2-4.0*a_ze*c_ze;
if Delta_ze==0
    ze=-b_ze/(2.0*a_ze);
elseif Delta_ze>0.0
    ze1=(-b_ze+sqrt(Delta_ze))/(2.0*a_ze);
    ze2=(-b_ze-sqrt(Delta_ze))/(2.0*a_ze);
    if ze1>=0.0 && ze2<0.0
        ze=ze1;
    elseif ze1<0.0 && ze2>=0.0
        ze=ze2;
    else
        disp('Unrealistic: roots of quadratic function of ze are both negative or positive, check it!')
    end
else
    disp('Unrealistic: delta of quadratic function is negative, no real value for ze, check it!')
end
end
%===function #2 for elastic solution===
function [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=SCeClayfmYtMt(zeTail,fm,n,R,Q,h,Lco,kint,yt,Mt)
E1e=R;
E2e=Q*h*h-4*R;
E3e=6*R-2*Q*h*h;
itrmaxe=1000;
tole=1.0e-5;
xpt=zeros(n+5,1);
x=zeros(n+5,1);
ye=zeros(n+5,1);
pe=zeros(n+5,1);
Ke=zeros(n+5,1);
Me=zeros(n+5,1);
Ve=zeros(n+5,1);
Se=zeros(n+5,1);
Fe=zeros(n+5,1);
Je=zeros(n+5,n+5);
for i=3:n+3
    xpt(i)=h*(i-3);
    x(i)=xpt(i)-Lco;
end
for itre=1:itrmaxe
    for i=3:n+3
        if zeTail(i)>0 
           pe(i)=fm*kint*ye(i);
           Ke(i)=fm*kint;
        elseif zeTail(i)==0
            pe(i)=0;
            Ke(i)=0;
        end 
        Fe(i)=E1e*ye(i+2)+E2e*ye(i+1)+E3e*ye(i)+h^4*pe(i)+E2e*ye(i-1)+E1e*ye(i-2);
        Je(i,i+2)=E1e;
        Je(i,i+1)=E2e;
        Je(i,i-1)=E2e;
        Je(i,i-2)=E1e;
        Je(i,i)=E3e+h^4*Ke(i);
    end 
    Fe(1)=ye(3)-yt;
    Je(1,3)=1;
    Fe(2)=R/h/h*(ye(4)-2*ye(3)+ye(2))-Mt;
    Je(2,4)=R/h/h;
    Je(2,3)=-2*R/h/h;
    Je(2,2)=R/h/h;
    Fe(n+4)=ye(n+4)-2*ye(n+3)+ye(n+2);
    Je(n+4,n+4)=1;
    Je(n+4,n+3)=-2;
    Je(n+4,n+2)=1;
    Fe(n+5)=ye(n+5)-2*ye(n+4)+2*ye(n+2)-ye(n+1);  
    Je(n+5,n+5)=1;
    Je(n+5,n+4)=-2;
    Je(n+5,n+2)=2;
    Je(n+5,n+1)=-1;
    deltae=-Je\Fe;
    ye=ye+deltae;
    maxdlte=max(deltae);
    erre=norm(deltae);
    if maxdlte<tole 
        disp(['[ELASTIC]Present loading step: iteration=',num2str(itre),' times, maximum delta=',num2str(maxdlte),', Euclidean norm of error=',num2str(erre)])
        break
    elseif maxdlte>tole && itre==itrmaxe
        warning('Iteration Failed...')
    end
end 
for i=3:n+3
    Me(i)=R/(h^2)*(ye(i+1)-2*ye(i)+ye(i-1));%(kN-m)
    Ve(i)=R/2/(h^3)*(-ye(i-2)+2*ye(i-1)-2*ye(i+1)+ye(i+2));%(kN)
    Se(i)=(ye(i+1)-ye(i-1))/2/h;%(rad)
end
outputx=x(3:n+3);%depth below pre-scour ground surface(m)
outputye_notailor=ye;%deflection from i=1 to i=n+5(including four imaginary nodes)
outputye=ye(3:n+3);%deflection(m)(i=3:n+3 without four imaginary nodes)
outputMe=Me(3:n+3);%bending moment(kN-m)
outputVe=Ve(3:n+3);%shear force(kN)
outputSe=rad2deg(Se(3:n+3));%rotation(degree)
outputpe=pe(3:n+3);%soil resistance(kN/m)
end
%===function #3 for Matlock solution===
function [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=SCMatlockfmYtMt(zeTail,fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt,Mt,ye)
%inputs:
%zeTail=n+5-by-1 vector saving equivalent depths considering scour-hole dim
%fm=p-multiplier for considering pile group effect
%D=pile diameter(m)
%n=number of pile elements
%R=flexural stiffness of pile, R=EI(kN-m2)
%Q=axial load at pile top(kN)
%h=element size(m)
%Lco=pre-scour cantilever length of pile(m)
%uw=buoyant unit weight of soil(kN/m3) as water table is above ground line
%cu=undrained shear strength of soft clay(kPa)
%Jclay=0.5 for soft clay (Matlock p-y curve)
%y50 is deflection (m) corresponding to p=pu/2 on Matlock p-y curve
%yt=deformation boundary condition: deflection applied at pile top(m)
%Mt=force boundary condition: bending moment applied at pile top(kN-m)
%ye=elastic solution of y(m), which is used as an initial guess in Matlock's iteration
xpt=zeros(n+5,1);
x=zeros(n+5,1);
for i=3:n+3
    xpt(i)=h*(i-3);%depth below pile head(m)
    x(i)=xpt(i)-Lco;%depth below pre-scour ground(m)
end
pu=zeros(n+5,1);
pus=zeros(n+5,1);
pud=zeros(n+5,1);
p=zeros(n+5,1);
F=zeros(n+5,1);
J=zeros(n+5,n+5);
k=zeros(n+5,1);
K=zeros(n+5,1);
E1=R;
E2=Q*h*h-4*R;
E3=6*R-2*Q*h*h;
itrmax=1000;
tol=1.0e-3;
y=ye;%use elastic response as the initial guess for y(m)
for i=3:n+3
    if zeTail(i)>0 
        pus(i)=(3+uw*zeTail(i)/cu+Jclay*zeTail(i)/D)*cu*D;%(kN/m)
        pud(i)=9*cu*D;%(kN/m)
        pu(i)=min(pus(i),pud(i));%(kN/m)
    else 
    end
end
for itr=1:itrmax
    for i=3:n+3
        if zeTail(i)>0 
            if y(i)<8*y50        
                p(i)=fm*0.5*pu(i)*nthroot(y(i)/y50,3);
                K(i)=fm*1/6*pu(i)*nthroot(1/y50,3)*nthroot(1/y(i)^2,3);
            else 
                p(i)=fm*pu(i);
                K(i)=0;
            end
        elseif zeTail(i)==0 
            p(i)=0;
            K(i)=0;
        end 
        F(i)=E1*y(i+2)+E2*y(i+1)+E3*y(i)+h^4*p(i)+E2*y(i-1)+E1*y(i-2);
        J(i,i+2)=E1;
        J(i,i+1)=E2;
        J(i,i-1)=E2;
        J(i,i-2)=E1;
        J(i,i)=E3+h^4*K(i);
    end 
    F(1)=y(3)-yt;
    J(1,3)=1;
    F(2)=R/h/h*(y(4)-2*y(3)+y(2))-Mt;
    J(2,4)=R/h/h;
    J(2,3)=-2*R/h/h;
    J(2,2)=R/h/h;
    F(n+4)=y(n+4)-2*y(n+3)+y(n+2);
    J(n+4,n+4)=1;
    J(n+4,n+3)=-2;
    J(n+4,n+2)=1;
    F(n+5)=y(n+5)-2*y(n+4)+2*y(n+2)-y(n+1);  
    J(n+5,n+5)=1;
    J(n+5,n+4)=-2;
    J(n+5,n+2)=2;
    J(n+5,n+1)=-1;
    delta=-J\F;
    y=y+delta;
    maxdlt=max(delta);
    err=norm(delta);
    if maxdlt<tol 
        disp(['[Matlock]Present loading step: iteration=',num2str(itr),' times, maximum delta=',num2str(maxdlt), ', Euclidean norm of error=',num2str(err)])
        break
    elseif maxdlt>tol && itr==itrmax
        warning('Iteration Failed...')
    end
end 
M=zeros(n+5,1);%bending moment
V=zeros(n+5,1);%shear force
S=zeros(n+5,1);%rotation
for i=3:n+3
    M(i)=R/(h^2)*(y(i+1)-2*y(i)+y(i-1));%(kN-m)
    V(i)=R/2/(h^3)*(-y(i-2)+2*y(i-1)-2*y(i+1)+y(i+2));%(kN)
    S(i)=(y(i+1)-y(i-1))/2/h;%(rad)
    k(i)=p(i)/y(i);%k=p/y(kPa)
end
outputy=y(3:n+3);%deflection from pile top to tip(m)
outputM=M(3:n+3);%bending moment(kN-m)
outputV=V(3:n+3);%shear force(kN)
outputS=rad2deg(S(3:n+3));%rotation(degree)
outputp=p(3:n+3);%soil resistance(kN/m)
outputpu=fm.*pu(3:n+3);%ultimate soil resistance(kN/m)
outputk=k(3:n+3);%k=p/y(kPa)
end
%==============================IN THE END====================================