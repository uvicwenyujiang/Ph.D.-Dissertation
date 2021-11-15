%*************************************************************************%
% This script was developed by Wenyu Jiang and Cheng Lin                  %
% Department of Civil Engineering                                         %
% University of Victoria, Canada                                          %
% Email: wenyujiang@uvic.ca; chenglin918@uvic.ca                          %
% Date: 2021-06-26                                                        %
%*************************************************************************%
%*************************************************************************%
% The code is for analyzing an individual pile in a group under laterally %
% loaded conditions. Soil-pile interactions are simulated using Matlock's %
% p-y curve.                                                              %  
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
%pile parameters:
R=28600.0;%flexural rigidity,R=EI(kN-m2)
L=12.07;%total length(m)
Lco=0.495;%unsupported length(m)
Leo=L-Lco;%embedded length(m)
D=0.324;%pile diameter(m)
Q=0;%axial load applied at pile head(kN)
fm=0.6;%p-multiplier to account for geotechnical pile group effect
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
%complete defining input values for computation.
%==================================================================================================
yeNOtailor_allstep=zeros(n+5,nstep);%(elastic)lateral displacement including four imaginary nodes(m)
ye_allstep=zeros(n+1,nstep);%(elastic)lateral displacement(m)
Me_allstep=zeros(n+1,nstep);%(elastic)bending moment(kN-m)
Ve_allstep=zeros(n+1,nstep);%(elastic)shear force(kN)
Se_allstep=zeros(n+1,nstep);%(elastic)rotation(deg)
pe_allstep=zeros(n+1,nstep);%(elastic)soil resistance(kN/m)
for loadstep=1:nstep
    [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=...
        eClayfmYtMt(fm,n,R,Q,h,Lco,kint,yt_vector(loadstep),Mt_vector(loadstep));
    yeNOtailor_allstep(:,loadstep)=outputye_notailor;
    ye_allstep(:,loadstep)=outputye;
    Me_allstep(:,loadstep)=outputMe;
    Ve_allstep(:,loadstep)=outputVe;
    Se_allstep(:,loadstep)=outputSe;
    pe_allstep(:,loadstep)=outputpe;
end %complete calculating elastic solutions
xnorm=outputx./D;%normalized depth below ground(xD)
y_allstep=zeros(n+1,nstep);%lateral displacement(m)
M_allstep=zeros(n+1,nstep);%bending moment(kN-m)
V_allstep=zeros(n+1,nstep);%shear force(kN)
S_allstep=zeros(n+1,nstep);%rotation(deg)
p_allstep=zeros(n+1,nstep);%soil resistance(kN/m)
k_allstep=zeros(n+1,nstep);%k=p/y (kPa)
pu_allstep=zeros(n+1);%ultimate soil resistance(kN/m)
for loadstep=1:nstep
    [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=...
        MatlockfmYtMt(fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt_vector(loadstep),Mt_vector(loadstep),yeNOtailor_allstep(:,loadstep));
    y_allstep(:,loadstep)=outputy;
    M_allstep(:,loadstep)=outputM;
    V_allstep(:,loadstep)=outputV;
    S_allstep(:,loadstep)=outputS;
    p_allstep(:,loadstep)=outputp;
    pu_allstep(:,1)=outputpu;
    k_allstep(:,loadstep)=outputk;
end
PTOPdisp=[0,y_allstep(1,:)]';%deflection at pile head(m)
PTOPshear=[0,V_allstep(1,:)]';%shear force at pile head(kN)
toc
%plot the calculated results:
figure(1)
subplot(2,4,1)
box on
grid on
plot(y_allstep(:,end),xnorm,'-r')
axis ij
xlabel('y(m)')
ylabel('Depth below ground(xD)')
title('Deflection')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,2)
box on
grid on
plot(S_allstep(:,end),xnorm,'-r')
axis ij
xlabel('\theta(^{\circ})')
ylabel('Depth below ground(xD)')
title('Rotation')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,3)
box on
grid on
plot(M_allstep(:,end),xnorm,'-r')
axis ij
xlabel('M(kN-m)')
ylabel('Depth below ground(xD)')
title('Bending Moment')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,4)
box on
grid on
plot(V_allstep(:,end),xnorm,'-r')
axis ij
xlabel('V(kN)')
ylabel('Depth below ground')
title('Shear force')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,5)
box on
grid on
plot(p_allstep(:,end),xnorm,'-r')
axis ij
hold on
plot(pu_allstep(:,1),xnorm,'-b')
legend('p','p_{u}')
xlabel('p(kN/m)')
ylabel('Depth below ground(xD)')
title('Soil Resistance')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,6)
box on
grid on
plot(p_allstep(:,end)./y_allstep(:,end),xnorm,'-r')
axis ij
xlabel('k=p/y(kPa)')
ylabel('Depth below ground(xD)')
title('k=p/y')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,7)
box on
grid on
plot(PTOPdisp,PTOPshear,'-o')
xlabel('y(m)')
ylabel('V(kN)')
title('Lateral load-deflection curve at pile head')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,8)
box on
grid on
plot([0,abs(M_allstep(1,:))./R]',[0,abs(M_allstep(1,:))]','-o')
xlabel('Curvature')
ylabel('Bending moment(kN-m)')
title('Moment-curvature curve at pile head')
set(gca,'FontName','Times','FontSize',11)
%===function #1 for elastic solution===
function [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=eClayfmYtMt(fm,n,R,Q,h,Lco,kint,yt,Mt)
E1e=R;
E2e=Q*h*h-4*R;
E3e=6*R-2*Q*h*h;
itrmaxe=1000;%predefined maximum number of iterations
tole=1.0e-5;%predefined tolerance for convergence control
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
        if x(i)>=0 
           pe(i)=fm*kint*ye(i);
           Ke(i)=fm*kint;
        else 
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
    Me(i)=R/(h^2)*(ye(i+1)-2*ye(i)+ye(i-1));
    Ve(i)=R/2/(h^3)*(-ye(i-2)+2*ye(i-1)-2*ye(i+1)+ye(i+2));
    Se(i)=(ye(i+1)-ye(i-1))/2/h;
end
ye(1)=ye(3);
ye(2)=ye(3);
ye(n+4)=ye(n+3);
ye(n+5)=ye(n+3);
outputx=x(3:n+3);%depth below ground surface(m)
outputye_notailor=ye;%deflection from i=1 to i=n+5(including four imaginary nodes)
outputye=ye(3:n+3);%deflection(m)(i=3:n+3 without four imaginary nodes)
outputMe=Me(3:n+3);%bending moment(kN-m)
outputVe=Ve(3:n+3);%shear force(kN)
outputSe=rad2deg(Se(3:n+3));%rotation(deg)
outputpe=pe(3:n+3);%soil resistance(kN/m)
end 
%===function #2 for Matlock solution===
function [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=MatlockfmYtMt(fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt,Mt,ye)
%fm=p-multiplier for considering pile group effect
%D=pile diameter(m)
%n=number of pile elements
%R=flexural stiffness of pile, R=EI(kN-m2)
%Q=axial load at pile top(kN)
%h=element size(m)
%Lco=pre-scour cantilever length of pile(m)
%uw=submerged unit weight of soil(kN/m3) as water table is above ground line
%cu=undrained shear strength of soft clay(kPa)
%Jclay=0.5 for soft clay (Matlock p-y curve)
%y50 is deflection (m) corresponding to p=pu/2 on Matlock p-y curve
%Vt=force boundary condition: shear force applied at pile top(kN)
%Mt=force boundary condition: bending moment applied at pile top(kN-m)
%ye=elastic solution of y(m), which is used as an initial guess in Matlock's iteration
xpt=zeros(n+5,1);
x=zeros(n+5,1);
for i=3:n+3
    xpt(i)=h*(i-3);%depth below pile head(m)
    x(i)=xpt(i)-Lco;%depth below pre-scour ground(m)
end
pu=zeros(n+5,1);%ultimate soil resistance of Matlock p-y curve(kN/m)
pus=zeros(n+5,1);%ultimate soil resistance of Matlock p-y curve(kN/m):shallow mode
pud=zeros(n+5,1);%ultimate soil resistance of Matlock p-y curve(kN/m):deep mode
p=zeros(n+5,1);%soil resistance dependent on y from p-y relation(kN/m)
F=zeros(n+5,1);%functions
J=zeros(n+5,n+5);%Jacobian matrix
k=zeros(n+5,1);%k=p/y not K=dp/dy
K=zeros(n+5,1);%K=dp/dy, tangent of a node on p-y curve
%define constants
E1=R;
E2=Q*h*h-4*R;
E3=6*R-2*Q*h*h;
itrmax=1000;%max. number of iterations
tol=1.0e-3;%predefined tolerance for convergence control
y=ye;%use elastic response as an initial guess for y:
for i=3:n+3
    if x(i)>=0 
        pus(i)=(3+uw*x(i)/cu+Jclay*x(i)/D)*cu*D;
        pud(i)=9*cu*D;
        pu(i)=min(pus(i),pud(i));
    else 
    end
end
for itr=1:itrmax 
    for i=3:n+3
        if x(i)>=0 
            if y(i)<8*y50        
                p(i)=fm*0.5*pu(i)*nthroot(y(i)/y50,3);
                K(i)=fm*1/6*pu(i)*nthroot(1/y50,3)*nthroot(1/y(i)^2,3);
            else 
                p(i)=fm*pu(i);
                K(i)=0;
            end       
        else 
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
outputS=rad2deg(S(3:n+3));%rotation (degree)
outputp=p(3:n+3);%soil resistance(kN/m)
outputpu=fm.*pu(3:n+3);%ultimate soil resistance(kN/m)
outputk=k(3:n+3);%k=p/y(kPa)
end
%=============IN THE END===================================================