%*************************************************************************%
% This script was developed by Wenyu Jiang and Cheng Lin                  %
% Department of Civil Engineering                                         %
% University of Victoria, Canada                                          %
% Email: wenyujiang@uvic.ca; chenglin918@uvic.ca                          %
% Date: 2021-06-26                                                        %
%*************************************************************************%
%*************************************************************************%
% The code is for analyzing a trailing row pile in a group under local    %
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
Ncol=2;%no. of piles on x-dir.
Nrow=2;%no. of piles on y-dir. 
D=0.324;%pile diameter(m)
ax=3.92*D;%pile spacing on x-dir.(m)
ay=3.29*D;%pile spacing on y-dir.(m)
R=28600.0;%flexural rigidity of pile,R=EI(kN-m2)
L=12.07;%pile length(m)
Lco=0.495;%unsupported length(m)
Leo=L-Lco;%embedded length(m)
Q=0;%axial load at pile head(kN)
prompt = 'What is row no. of the pile in the group? ';%x-dir
prowinput = input(prompt)
prompt = 'What is col. no. of the pile in the group? ';%y-dir
pcolinput = input(prompt)
fm=0.6;%p-multiplier to account for geotechnical group effect
SSIeta=1.0;%correction factor to account for soil-pile interaction
%scour-hole parameters
Sd=1.5*D;%6.0593*D;%scour-hole depth(m)
Swb=0*D;%scour-hole bottom width(m)
Sq=26.6;%scour-hole slope angle(deg)
Swt=Swb+Sd*cotd(Sq);%scour-hole top width(m)
scbtm_x=(Ncol-1)*ax+2*Swb+D;%length of scour-hole bottom(m)
scbtm_y=(Nrow-1)*ay+2*Swb+D;%width of scour-hole bottom(m)
scbtm_area_sym=scbtm_x*scbtm_y-4*((Swb+D/2)^2)+sym(pi)*((Swb+D/2)^2);%area of scour-hole bottom(m2)
scbtm_area=double(scbtm_area_sym);%area of scour-hole bottom(m2)
%Swc=distance between circumference of center pile and the periphery of the circlular scour-hole
%bottom (used in the closed-form analytical method)
Swc=double(sqrt(scbtm_area_sym/sym(pi))-D/2);%(m)
lwratio=scbtm_x/scbtm_y;%length-to-width ratio of scour bottom
Ve=sqrt(scbtm_area/lwratio);%width(y) equivalent rectanglar scour bottom(m)
Ue=lwratio*Ve;%legnth(x) equivalent rectanglar scour bottom(m)
%post-scour pile unsupported and embedded lengths
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
%--------------------------------------------------------------------------
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
        zblwsc(page)=0.0001;%(m)
    else
    end
end
zblwscnorm=zblwsc./D;%normalized depth below scour bottom(xD)
%mesh generation for calculating the equivalent depth
dx1dsr=0.015*D;%desired size of delta x1(m)
n1=ceil(Sd*cotd(Sq)/dx1dsr);%element number n1
dx2dsr=0.015*D;%desired size of delta x2(m)
n2=ceil(Ue/dx2dsr);%element number n2
dx3dsr=0.015*D;%desired size of delta x3(m)
n3=ceil(Ve/dx3dsr);%element number n3
%middle-middle zone is for piles
xpndbox=zeros(Nrow,Ncol);
ypndbox=zeros(Nrow,Ncol);
%upper-left pile is the reference pile
xpref=D/2+Swb+(Ue-scbtm_x)/2;%(m)
ypref=D/2+Swb+(Ve-scbtm_y)/2;%(m)
for row=1:Nrow
    xpndbox(row,:)=linspace(xpref,xpref+(Ncol-1)*ax,Ncol);%(m)
    ypndbox(row,:)=linspace(ypref+(row-1)*ay,ypref+(row-1)*ay,Ncol);%(m)
end
xpnd=xpndbox(prowinput,pcolinput);
ypnd=ypndbox(prowinput,pcolinput);
%interior slope sidewall regions(zone1-zone8)
%1.upper-left corner zone
nx1=n1;
ny1=n1;
dx1=Sd*cotd(Sq)/nx1;
dy1=Sd*cotd(Sq)/ny1;
xf1=zeros(ny1,nx1);
yf1=zeros(ny1,nx1);
h1=zeros(ny1,nx1);
dA1=zeros(ny1,nx1);
for row=1:ny1
    for col=1:nx1
        if col<=row
            h1(row,col)=(2*row-1)*dx1/2*tand(Sq);
        else
            h1(row,col)=(2*col-1)*dx1/2*tand(Sq);
        end
    end
end
for col=1:nx1
    xf1(:,col)=linspace((1-2*col)*dx1/2,(1-2*col)*dx1/2,ny1)';%(m)
    yf1(:,col)=linspace(-dy1/2,(1-2*ny1)*dy1/2,ny1)';%(m)
    dA1(:,col)=linspace(dx1*dy1,dx1*dy1,ny1)';%(m2)
end
dv1=dA1.*h1;%(m3)
%2.upper-middle zone
nx2=n2;
ny2=n1;
dx2=Ue/nx2;
dy2=Sd*cotd(Sq)/ny2;
xf2=zeros(ny2,nx2);
yf2=zeros(ny2,nx2);
h2=zeros(ny2,nx2);
dA2=zeros(ny2,nx2);
for row=1:ny2
    xf2(row,:)=linspace(dx2/2,(2*nx2-1)*dx2/2,nx2);
    yf2(row,:)=linspace((1-2*row)*dy2/2,(1-2*row)*dy2/2,nx2);
    dA2(row,:)=linspace(dx2*dy2,dx2*dy2,nx2);
    h2(row,:)=linspace((2*row-1)*dy2/2*tand(Sq),(2*row-1)*dy2/2*tand(Sq),nx2);
    dv2(row,:)=linspace(dx2*dy2*(2*row-1)*dy2/2*tand(Sq),...
        dx2*dy2*(2*row-1)*dy2/2*tand(Sq),nx2);
end
%3.upper-right corner zone
nx3=n1;
ny3=n1;
dx3=Sd*cotd(Sq)/nx3;
dy3=Sd*cotd(Sq)/ny3;
xf3=zeros(ny3,nx3);
yf3=zeros(ny3,nx3);
h3=zeros(ny3,nx3);
dA3=zeros(ny3,nx3);
for row=1:ny3
    for col=1:nx3
        if col<=row
            h3(row,col)=(2*row-1)*dx3/2*tand(Sq);
        else
            h3(row,col)=(2*col-1)*dx3/2*tand(Sq);
        end
    end
end
for col=1:nx3
    xf3(:,col)=linspace(Ue+(2*col-1)*dx3/2,Ue+(2*col-1)*dx3/2,ny3)';%(m)
    yf3(:,col)=linspace(-dy3/2,(1-2*ny3)*dy3/2,ny3)';%(m)
    dA3(:,col)=linspace(dx3*dy3,dx3*dy3,ny3)';%(m2)
end
dv3=dA3.*h3;%(m3)
%4.middle-left zone
nx4=n1;
ny4=n3;
dx4=Sd*cotd(Sq)/nx4;
dy4=Ve/ny4;
xf4=zeros(ny4,nx4);
yf4=zeros(ny4,nx4);
h4=zeros(ny4,nx4);
dA4=zeros(ny4,nx4);
dv4=zeros(ny4,nx4);
for col=1:nx4
    xf4(:,col)=linspace((1-2*col)*dx4/2,(1-2*col)*dx4/2,ny4)';
    yf4(:,col)=linspace(dy4/2,(2*ny4-1)*dy4/2,ny4)';
    h4(:,col)=linspace((2*col-1)*dx4/2*tand(Sq),(2*col-1)*dx4/2*tand(Sq),ny4)';
    dA4(:,col)=linspace(dx4*dy4,dx4*dy4,ny4)';
    dv4(:,col)=linspace(dx4*dy4*(2*col-1)*dx4/2*tand(Sq),...
        dx4*dy4*(2*col-1)*dx4/2*tand(Sq),ny4)';
end
%5.middle-right zone
nx5=n1;
ny5=n3;
dxf5=Sd*cotd(Sq)/nx5;
dyf5=Ve/ny5;
xf5=zeros(ny5,nx5);
yf5=zeros(ny5,nx5);
h5=zeros(ny5,nx5);
dv5=zeros(ny5,nx5);
for col=1:nx5
    xf5(:,col)=linspace(Ue+(2*col-1)*dxf5/2,Ue+(2*col-1)*dxf5/2,ny5)';
    yf5(:,col)=linspace(dyf5/2,(2*ny5-1)*dyf5/2,ny5)';
    h5(:,col)=linspace((2*col-1)*dxf5/2*tand(Sq),(2*col-1)*dxf5/2*tand(Sq),ny5)';
    dA5(:,col)=linspace(dxf5*dyf5,dxf5*dyf5,ny5)';
    dv5(:,col)=linspace(dxf5*dyf5*(2*col-1)*dxf5/2*tand(Sq),...
        dxf5*dyf5*(2*col-1)*dxf5/2*tand(Sq),ny5)';
end
%6.lower-left corner zone
nx6=n1;
ny6=n1;
dx6=Sd*cotd(Sq)/nx6;
dy6=Sd*cotd(Sq)/ny6;
xf6=zeros(ny6,nx6);
yf6=zeros(ny6,nx6);
h6=zeros(ny6,nx6);
dA6=zeros(ny6,nx6);
for row=1:ny6
    for col=1:nx6
        if col<=row
            h6(row,col)=(2*row-1)*dx6/2*tand(Sq);
        else
            h6(row,col)=(2*col-1)*dx6/2*tand(Sq);
        end
    end
end
for col=1:nx6
    xf6(:,col)=linspace((1-2*col)*dx6/2,(1-2*col)*dx6/2,ny6)';%(m)
    yf6(:,col)=linspace(Ve+dy6/2,Ve+(2*ny6-1)*dy6/2,ny6)';%(m)
    dA6(:,col)=linspace(dx6*dy6,dx6*dy6,ny6)';%(m2)
end
dv6=dA6.*h6;%(m3)
%7.lower-middle zone
nx7=n2;
ny7=n1;
dx7=Ue/nx7;
dy7=Sd*cotd(Sq)/ny7;
xf7=zeros(ny7,nx7);
yf7=zeros(ny7,nx7);
h7=zeros(ny7,nx7);
dA7=zeros(ny7,nx7);
dv7=zeros(ny7,nx7);
for row=1:ny7
    xf7(row,:)=linspace(dx7/2,(2*nx7-1)*dx7/2,nx7);
    yf7(row,:)=linspace(Ve+(2*row-1)*dy7/2,Ve+(2*row-1)*dy7/2,nx7);
    dA7(row,:)=linspace(dx7*dy7,dx7*dy7,nx7);
    h7(row,:)=linspace((2*row-1)*dy7/2*tand(Sq),(2*row-1)*dy7/2*tand(Sq),nx7);
    dv7(row,:)=linspace(dx7*dy7*(2*row-1)*dy7/2*tand(Sq),...
        dx7*dy7*(2*row-1)*dy7/2*tand(Sq),nx7);
end
%8.lower-right corner zone
nx8=n1;
ny8=n1;
dx8=Sd*cotd(Sq)/nx8;
dy8=Sd*cotd(Sq)/ny8;
xf8=zeros(ny8,nx8);
yf8=zeros(ny8,nx8);
h8=zeros(ny8,nx8);
dA8=zeros(ny8,nx8);
for row=1:ny8
    for col=1:nx8
        if col<=row
            h8(row,col)=(2*row-1)*dx8/2*tand(Sq);
        else
            h8(row,col)=(2*col-1)*dx8/2*tand(Sq);
        end
    end
end
for col=1:nx8
    xf8(:,col)=linspace(Ue+(2*col-1)*dx8/2,Ue+(2*col-1)*dx8/2,ny8)';%(m)
    yf8(:,col)=linspace(Ve+dy8/2,Ve+(2*ny8-1)*dy8/2,ny8)';%(m)
    dA8(:,col)=linspace(dx8*dy8,dx8*dy8,ny8)';%(m2)
end
dv8=dA8.*h8;%(m3)
%compelete defining interior sidewall loading elements, and mesh for exterior regions(zone1-zone8)
ex=200*D;
ey=200*D;
%element size and number for zone E1
dxE1dsr=0.5*D;%desired size of delta xE1(m)
nxE1=ceil((ex+Sd*cotd(Sq))/dxE1dsr);%element number nE1
dyE1dsr=0.5*D;%desired size of delta xE1(m)
nyE1=ceil(ey/dyE1dsr);%element number nE1
%element size and number for zone E2
dxE2dsr=0.5*D;%desired size of delta xE2(m)
nxE2=ceil(Ue/dxE2dsr);%element number nE2
dyE2dsr=0.5*D;%desired size of delta xE2(m)
nyE2=ceil(ey/dyE2dsr);%element number nE2
%element size and number for zone E3
dxE3dsr=0.5*D;%desired size of delta xE3(m)
nxE3=ceil((ex+Sd*cotd(Sq))/dxE3dsr);%element number nE3
dyE3dsr=0.5*D;%desired size of delta xE3(m)
nyE3=ceil(ey/dyE3dsr);%element number nE3
%element size and number for zone E4
dxE4dsr=0.5*D;%desired size of delta xE4(m)
nxE4=ceil(ex/dxE4dsr);%element number nE4
dyE4dsr=0.5*D;%desired size of delta xE4(m)
nyE4=ceil((Ve+2*Sd*cotd(Sq))/dyE4dsr);%element number nE4
%element size and number for zone E5
dxE5dsr=0.5*D;%desired size of delta xE5(m)
nxE5=ceil(ex/dxE5dsr);%element number nE5
dyE5dsr=0.5*D;%desired size of delta xE5(m)
nyE5=ceil((Ve+2*Sd*cotd(Sq))/dyE5dsr);%element number nE5
%element size and number for zone E6
dxE6dsr=0.5*D;%desired size of delta xE6(m)
nxE6=ceil((ex+Sd*cotd(Sq))/dxE6dsr);%element number nE6
dyE6dsr=0.5*D;%desired size of delta xE6(m)
nyE6=ceil(ey/dyE6dsr);%element number nE6
%element size and number for zone E7
dxE7dsr=0.5*D;%desired size of delta xE7(m)
nxE7=ceil(Ue/dxE7dsr);%element number nE7
dyE7dsr=0.5*D;%desired size of delta xE7(m)
nyE7=ceil(ey/dyE7dsr);%element number nE7
%element size and number for zone E8
dxE8dsr=0.5*D;%desired size of delta xE8(m)
nxE8=ceil((ex+Sd*cotd(Sq))/dxE8dsr);%element number nE8
dyE8dsr=0.5*D;%desired size of delta xE8(m)
nyE8=ceil(ey/dyE8dsr);%element number nE8
%1.upper-left zone E1
dxE1=(ex+Sd*cotd(Sq))/nxE1;
dyE1=ey/nyE1;
xE1=zeros(nyE1,nxE1);
yE1=zeros(nyE1,nxE1);
for row=1:nyE1
    xE1(row,:)=linspace(-dxE1/2,(1-2*nxE1)*dxE1/2,nxE1);
    yE1(row,:)=linspace((1-2*row)*dyE1/2-Sd*cotd(Sq),...
        (1-2*row)*dyE1/2-Sd*cotd(Sq),nxE1);
end
hE1=Sd.*ones(nyE1,nxE1);
dAE1=(dxE1*dyE1).*ones(nyE1,nxE1);
dvE1=dAE1.*hE1;
%2.upper-middle zone E2
dxE2=Ue/nxE2;
dyE2=ey/nyE2;
xE2=zeros(nyE2,nxE2);
yE2=zeros(nyE2,nxE2);
for row=1:nyE2
    xE2(row,:)=linspace(dxE2/2,(2*nxE2-1)*dxE2/2,nxE2);
    yE2(row,:)=linspace((1-2*row)*dyE2/2-Sd*cotd(Sq),...
        (1-2*row)*dyE2/2-Sd*cotd(Sq),nxE2);
end
hE2=Sd.*ones(nyE2,nxE2);
dAE2=(dxE2*dyE2).*ones(nyE2,nxE2);
dvE2=dAE2.*hE2;
%3.upper-right zone E3
dxE3=(ex+Sd*cotd(Sq))/nxE3;
dyE3=ey/nyE3;
xE3=zeros(nyE3,nxE3);
yE3=zeros(nyE3,nxE3);
for row=1:nyE3
    xE3(row,:)=linspace(Ue+dxE3/2,Ue+(2*nxE3-1)*dxE3/2,nxE3);
    yE3(row,:)=linspace((1-2*row)*dyE3/2-Sd*cotd(Sq),...
        (1-2*row)*dyE3/2-Sd*cotd(Sq),nxE3);
end
hE3=Sd.*ones(nyE3,nxE3);
dAE3=(dxE3*dyE3).*ones(nyE3,nxE3);
dvE3=dAE3.*hE3;
%4.middle-left zone E4
dxE4=ex/nxE4;
dyE4=(2*Sd*cotd(Sq)+Ve)/nyE4;
xE4=zeros(nyE4,nxE4);
yE4=zeros(nyE4,nxE4);
for col=1:nxE4
    xE4(:,col)=linspace((1-2*col)*dxE4/2-Sd*cotd(Sq),...
        (1-2*col)*dxE4/2-Sd*cotd(Sq),nyE4)';
    yE4(:,col)=linspace(dyE4/2-Sd*cotd(Sq),Ve+Sd*cotd(Sq)-dyE4/2,nyE4)';
end
hE4=Sd.*ones(nyE4,nxE4);
dAE4=(dxE4*dyE4).*ones(nyE4,nxE4);
dvE4=dAE4.*hE4;
%5.middle-right zone E5
dxE5=ex/nxE5;
dyE5=(2*Sd*cotd(Sq)+Ve)/nyE5;
xE5=zeros(nyE5,nxE5);
yE5=zeros(nyE5,nxE5);
for col=1:nxE5
    xE5(:,col)=linspace(Ue+Sd*cotd(Sq)+(2*col-1)*dxE5/2,...
        Ue+Sd*cotd(Sq)+(2*col-1)*dxE5/2,nyE5)';
    yE5(:,col)=linspace(dyE5/2-Sd*cotd(Sq),Ve+Sd*cotd(Sq)-dyE5/2,nyE5)';
end
hE5=Sd.*ones(nyE5,nxE5);
dAE5=(dxE5*dyE5).*ones(nyE5,nxE5);
dvE5=dAE5.*hE5;
%6.lower-left zone E6
dxE6=(ex+Sd*cotd(Sq))/nxE6;
dyE6=ey/nyE6;
xE6=zeros(nyE6,nxE6);
yE6=zeros(nyE6,nxE6);
for row=1:nyE6
    xE6(row,:)=linspace(-dxE6/2,(1-2*nxE6)*dxE6/2,nxE6);
    yE6(row,:)=linspace(Ve+Sd*cotd(Sq)+(2*row-1)*dyE6/2,...
        Ve+Sd*cotd(Sq)+(2*row-1)*dyE6/2,nxE6);
end
hE6=Sd.*ones(nyE6,nxE6);
dAE6=(dxE6*dyE6).*ones(nyE6,nxE6);
dvE6=dAE6.*hE6;
%7.lower-middle zone E7
dxE7=Ue/nxE7;
dyE7=ey/nyE7;
xE7=zeros(nyE7,nxE7);
yE7=zeros(nyE7,nxE7);
for row=1:nyE7
    xE7(row,:)=linspace(dxE7/2,(2*nxE7-1)*dxE7/2,nxE7);
    yE7(row,:)=linspace(Ve+Sd*cotd(Sq)+(2*row-1)*dyE7/2,...
        Ve+Sd*cotd(Sq)+(2*row-1)*dyE7/2,nxE7);
end
hE7=Sd.*ones(nyE7,nxE7);
dAE7=(dxE7*dyE7).*ones(nyE7,nxE7);
dvE7=dAE7.*hE7;
%8.lower-right zone E8
dxE8=(ex+Sd*cotd(Sq))/nxE8;
dyE8=ey/nyE8;
xE8=zeros(nyE8,nxE8);
yE8=zeros(nyE8,nxE8);
for row=1:nyE8
    xE8(row,:)=linspace(Ue+dxE8/2,Ue+(2*nxE8-1)*dxE8/2,nxE8);
    yE8(row,:)=linspace(Ve+Sd*cotd(Sq)+(2*row-1)*dyE8/2,...
        Ve+Sd*cotd(Sq)+(2*row-1)*dyE8/2,nxE8);
end
hE8=Sd.*ones(nyE8,nxE8);
dAE8=(dxE8*dyE8).*ones(nyE8,nxE8);
dvE8=dAE8.*hE8;
%calculate point loads:
p1=uw.*dv1;%(kN)
dlt1=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz1=zeros(ny1,nx1);
    for row=1:ny1
        for col=1:nx1
            dltatz1(row,col)=3*p1(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf1(row,col))^2+(ypnd-yf1(row,col))^2))^5);%(kPa)
        end
    end
    dlt1(page,1)=sum(dltatz1,'all');
end
%2.upper middle zone
p2=uw.*dv2;%(kN)
dlt2=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz2=zeros(ny2,nx2);
    for row=1:ny2
        for col=1:nx2
            dltatz2(row,col)=3*p2(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf2(row,col))^2+(ypnd-yf2(row,col))^2))^5);%(kPa)
        end
    end
    dlt2(page,1)=sum(dltatz2,'all');
end
%3.upper right zone
p3=uw.*dv3;%(kN)
dlt3=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz3=zeros(ny3,nx3);
    for row=1:ny3
        for col=1:nx3
            dltatz3(row,col)=3*p3(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf3(row,col))^2+(ypnd-yf3(row,col))^2))^5);%(kPa)
        end
    end
    dlt3(page,1)=sum(dltatz3,'all');
end
%4.middle left zone
p4=uw.*dv4;%(kN)
dlt4=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz4=zeros(ny4,nx4);
    for row=1:ny4
        for col=1:nx4
            dltatz4(row,col)=3*p4(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf4(row,col))^2+(ypnd-yf4(row,col))^2))^5);%(kPa)
        end
    end
    dlt4(page,1)=sum(dltatz4,'all');
end
%5.middle right zone
p5=uw.*dv5;%(kN)
dlt5=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz5=zeros(ny5,nx5);
    for row=1:ny5
        for col=1:nx5
            dltatz5(row,col)=3*p5(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf5(row,col))^2+(ypnd-yf5(row,col))^2))^5);%(kPa)
        end
    end
    dlt5(page,1)=sum(dltatz5,'all');
end
%6.lower left zone
p6=uw.*dv6;%(kN)
dlt6=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz6=zeros(ny6,nx6);
    for row=1:ny6
        for col=1:nx6
            dltatz6(row,col)=3*p6(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf6(row,col))^2+(ypnd-yf6(row,col))^2))^5);%(kPa)
        end
    end
    dlt6(page,1)=sum(dltatz6,'all');
end
%7.lower middle zone
p7=uw.*dv7;%(kN)
dlt7=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz7=zeros(ny7,nx7);
    for row=1:ny7
        for col=1:nx7
            dltatz7(row,col)=3*p7(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf7(row,col))^2+(ypnd-yf7(row,col))^2))^5);%(kPa)
        end
    end
    dlt7(page,1)=sum(dltatz7,'all');
end
%8.lower right zone
p8=uw.*dv8;%(kN)
dlt8=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatz8=zeros(ny8,nx8);
    for row=1:ny8
        for col=1:nx8
            dltatz8(row,col)=3*p8(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xf8(row,col))^2+(ypnd-yf8(row,col))^2))^5);%(kPa)
        end
    end
    dlt8(page,1)=sum(dltatz8,'all');
end
%1.upper left zone E1
pE1=uw.*dvE1;%(kN)
dltE1=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE1=zeros(nyE1,nxE1);
    for row=1:nyE1
        for col=1:nxE1
            dltatzE1(row,col)=3*pE1(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE1(row,col))^2+(ypnd-yE1(row,col))^2))^5);%(kPa)
        end
    end
    dltE1(page,1)=sum(dltatzE1,'all');
end
%2.upper middle zone E2
pE2=uw.*dvE2;%(kN)
dltE2=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE2=zeros(nyE2,nxE2);
    for row=1:nyE2
        for col=1:nxE2
            dltatzE2(row,col)=3*pE2(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE2(row,col))^2+(ypnd-yE2(row,col))^2))^5);%(kPa)
        end
    end
    dltE2(page,1)=sum(dltatzE2,'all');
end
%3.upper right zone E3
pE3=uw.*dvE3;%(kN)
dltE3=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE3=zeros(nyE3,nxE3);
    for row=1:nyE3
        for col=1:nxE3
            dltatzE3(row,col)=3*pE3(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE3(row,col))^2+(ypnd-yE3(row,col))^2))^5);%(kPa)
        end
    end
    dltE3(page,1)=sum(dltatzE3,'all');
end
%4.middle left zone E4
pE4=uw.*dvE4;%(kN)
dltE4=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE4=zeros(nyE4,nxE4);
    for row=1:nyE4
        for col=1:nxE4
            dltatzE4(row,col)=3*pE4(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE4(row,col))^2+(ypnd-yE4(row,col))^2))^5);%(kPa)
        end
    end
    dltE4(page,1)=sum(dltatzE4,'all');
end
%5.middle right zone E5
pE5=uw.*dvE5;%(kN)
dltE5=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE5=zeros(nyE5,nxE5);
    for row=1:nyE5
        for col=1:nxE5
            dltatzE5(row,col)=3*pE5(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE5(row,col))^2+(ypnd-yE5(row,col))^2))^5);%(kPa)
        end
    end
    dltE5(page,1)=sum(dltatzE5,'all');
end
%6.lower left zone E6
pE6=uw.*dvE6;%(kN)
dltE6=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE6=zeros(nyE6,nxE6);
    for row=1:nyE6
        for col=1:nxE6
            dltatzE6(row,col)=3*pE6(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE6(row,col))^2+(ypnd-yE6(row,col))^2))^5);%(kPa)
        end
    end
    dltE6(page,1)=sum(dltatzE6,'all');
end
%7.lower middle zone E7
pE7=uw.*dvE7;%(kN)
dltE7=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE7=zeros(nyE7,nxE7);
    for row=1:nyE7
        for col=1:nxE7
            dltatzE7(row,col)=3*pE7(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE7(row,col))^2+(ypnd-yE7(row,col))^2))^5);%(kPa)
        end
    end
    dltE7(page,1)=sum(dltatzE7,'all');
end
%8.lower right zone E8
pE8=uw.*dvE8;%(kN)
dltE8=zeros(Nndblwsc,1);
for page=1:Nndblwsc
    dltatzE8=zeros(nyE8,nxE8);
    for row=1:nyE8
        for col=1:nxE8
            dltatzE8(row,col)=3*pE8(row,col)*zblwsc(page)^3/...
           (2*pi*(sqrt(zblwsc(page)^2+(xpnd-xE8(row,col))^2+(ypnd-yE8(row,col))^2))^5);%(kPa)
        end
    end
    dltE8(page,1)=sum(dltatzE8,'all');
end
%calculate additional vertical effective stress: 
dsigmaZ=dlt1+dlt2+dlt3+dlt4+dlt5+dlt6+dlt7+dlt8+...
    dltE1+dltE2+dltE3+dltE4+dltE5+dltE6+dltE7+dltE8;%(kPa)
sigmaZsimp=uw.*zblwsc;%(kPa)
sigmaZsc=SSIeta.*(sigmaZsimp+dsigmaZ);%(kPa)
sigmaZint=uw.*(zblwsc+Sd);%(kPa)
for page=1:Nndblwsc
    if sigmaZsc(page)>sigmaZint(page) 
        sigmaZsc(page)=sigmaZint(page);
    else
    end
end
ze=sigmaZsc./uw;%equivalent depth(m)
SSR=sigmaZsc./sigmaZint;%vertical stress ratio
%Besides, closed-form analytical method:
dltAM=zeros(Nndblwsc,1);%(kPa)
for page=1:Nndblwsc
    dltAM(page,1)=uw*zblwsc(page)*tand(Sq)*...
        ((Sd/tand(Sq)+Swc)/sqrt((Sd/tand(Sq)+Swc)^2+zblwsc(page)^2)-...
        Swc/sqrt(Swc^2+zblwsc(page)^2));%(kPa)
end
sigmava_AM=sigmaZsimp+dltAM;%(kPa)
ze_AM=sigmava_AM./uw;%(m)
SSR_AM=sigmava_AM./sigmaZint;
%tailor ze-vectors:
ie1=(n+3)-Nndblwsc+1;%node tag of the 1st embedded node in a n+5 vector
zeTail=[linspace(0,0,ie1-1)';ze;0;0];%tailored vector(length=n+5)
%Linear elastic solution:
yeNOtailor_allstep=zeros(n+5,nstep);%elastic lateral displacement including four imaginary nodes' results(m)
ye_allstep=zeros(n+1,nstep);%elastic lateral displacement(m)
Me_allstep=zeros(n+1,nstep);%elastic bending moment(kN-m)
Ve_allstep=zeros(n+1,nstep);%elastic shear force(kN)
Se_allstep=zeros(n+1,nstep);%elastic rotation(deg)
pe_allstep=zeros(n+1,nstep);%elastic soil resistance(kN/m)
for loadstep=1:nstep
    [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=SCeClayfmYtMt(zeTail,fm,n,R,Q,h,Lco,kint,yt_vector(loadstep),Mt_vector(loadstep));
    yeNOtailor_allstep(:,loadstep)=outputye_notailor;%(m)
    ye_allstep(:,loadstep)=outputye;%(m)
    Me_allstep(:,loadstep)=outputMe;%(kN-m)
    Ve_allstep(:,loadstep)=outputVe;%(kN)
    Se_allstep(:,loadstep)=outputSe;%(deg)
    pe_allstep(:,loadstep)=outputpe;%(kN/m)
end
xnorm=outputx./D;%normalized depth below pre-scour ground(xD)
%Non-linear solution (Matlock's p-y curve):
y_allstep=zeros(n+1,nstep);%lateral displacement(m)
M_allstep=zeros(n+1,nstep);%bending moment(kN-m)
V_allstep=zeros(n+1,nstep);%shear force(kN)
S_allstep=zeros(n+1,nstep);%rotation(deg)
p_allstep=zeros(n+1,nstep);%soil resistance(kN/m)
k_allstep=zeros(n+1,nstep);%k=p/y (kPa)
pu_allstep=zeros(n+1,1);%ultimate soil resistance(kN/m)
for loadstep=1:nstep
    [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=...
        SCMatlockfmYtMt(zeTail,fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt_vector(loadstep),Mt_vector(loadstep),yeNOtailor_allstep(:,loadstep));
    y_allstep(:,loadstep)=outputy;%(m)
    M_allstep(:,loadstep)=outputM;%(kN-m)
    V_allstep(:,loadstep)=outputV;%(kN)
    S_allstep(:,loadstep)=outputS;%(deg)
    p_allstep(:,loadstep)=outputp;%(kN/m)
    pu_allstep(:,1)=outputpu;%(kN/m)
    k_allstep(:,loadstep)=outputk;%(kPa)
end
PTOPdisp=[0,y_allstep(1,:)]';%pile-top deflection(m)
PTOPshear=[0,V_allstep(1,:)]';%pile-top load(kN)
toc 
%Plot the calculated results of pile
figure(1)
subplot(2,4,1)
box on
grid on
plot(y_allstep(:,end),xnorm,'-r')
axis ij
xlabel('y(m)')
ylabel('Depth below pre-scour ground(xD)')
title('Deflection')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,2)
box on
grid on
plot(S_allstep(:,end),xnorm,'-r')
axis ij
xlabel('\theta(^{\circ})')
ylabel('Depth below pre-scour ground(xD)')
title('Rotation')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,3)
box on
grid on
plot(M_allstep(:,end),xnorm,'-r')
axis ij
xlabel('M(kN-m)')
ylabel('Depth below pre-scour ground(xD)')
title('Bending Moment')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,4)
box on
grid on
plot(V_allstep(:,end),xnorm,'-r')
axis ij
xlabel('V(kN)')
ylabel('Depth below pre-scour ground(xD)')
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
ylabel('Depth below pre-scour ground(xD)')
title('Soil resistance')
set(gca,'FontName','Times','FontSize',11)
subplot(2,4,6)
box on
grid on
plot(p_allstep(:,end)./y_allstep(:,end),xnorm,'-r')
axis ij
xlabel('k=p/y(kPa)')
ylabel('Depth below pre-scour ground(xD)')
title('Stiffness k=p/y')
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
title('Bending moment-curvature curve at pile head')
set(gca,'FontName','Times','FontSize',11)
%===function #1 for elastic solution===
function [outputx,outputye_notailor,outputye,outputMe,outputVe,outputSe,outputpe]=SCeClayfmYtMt(zeTail,fm,n,R,Q,h,Lco,kint,yt,Mt)
%inputs: 
%zeTail=n+5-by-1 vector saving equivalent depths considering scour-hole dim
%fm=p-multiplier considering pile group effects
%D=pile diameter(m)
%n=number of pile elements
%R=flexural stiffness of pile, R=EI(kN-m2)
%Q=axial load at pile top(kN)
%h=element size(m)
%Lco=pre-scour cantilever length of pile(m)
%kint=k(kN/m3) of initial portion of p-y curve(kN/m3)p=Es*y=(k*x)*y. kint=12600 kPa from CIVE580 project #2 [OneNote 2020-08-23: I found kint greater than 10000 kPa is necessary for convergency]
%yt=defomration boundary condition: deflection applied at pile top(m)
%Mt=force boundary condition: bending moment applied at pile top(kN-m)
E1e=R;
E2e=Q*h*h-4*R;
E3e=6*R-2*Q*h*h;
itrmaxe=1000;%max. number of iterations
tole=1.0e-5;%tolerance for convergence control
xpt=zeros(n+5,1);%depth below pile top(m)
x=zeros(n+5,1);%depth below pre-scour ground(m)
ye=zeros(n+5,1);%pile deflection(m)
pe=zeros(n+5,1);%soil resistance(kN/m)
Ke=zeros(n+5,1);%K=dp/dy,tangent of a node on p-y curve
Me=zeros(n+5,1);%bending moment
Ve=zeros(n+5,1);%shear force
Se=zeros(n+5,1);%rotation
Fe=zeros(n+5,1);
Je=zeros(n+5,n+5);
for i=3:n+3
    xpt(i)=h*(i-3);%depth below pile head(m)
    x(i)=xpt(i)-Lco;%depth below pre-scour ground(m)
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
outputye=ye(3:n+3);%deflection from pile top to tip(i=3:n+3 without four imaginary nodes)(m)
outputMe=Me(3:n+3);%bending moment(kN-m)
outputVe=Ve(3:n+3);%shear force(kN)
outputSe=rad2deg(Se(3:n+3));%rotation(degree)
outputpe=pe(3:n+3);%soil resistance(kN/m)
end
%===function #2 for Matlock solution===
function [outputp,outputy,outputM,outputV,outputS,outputpu,outputk]=SCMatlockfmYtMt(zeTail,fm,D,n,R,Q,h,Lco,uw,cu,Jclay,y50,yt,Mt,ye)
%inputs:
%zeTail=n+5-by-1 vector saving equivalent depths considering scour-hole dim
%fm=p-multiplier for considering pile group effect
%D=pile diameter(m)
%n=number of pile elements
%R=flexural stiffness of pile, R=EI(kN-m2)
%Q=axial load at pile top(kN)
%h=element size(m)
%Lco=pre-scour unsupported length of pile(m)
%uw=submerged unit weight of soil(kN/m3)
%cu=undrained shear strength of soft clay(kPa)
%Jclay=0.5 for soft clay (Matlock p-y curve)
%y50 is deflection (m) corresponding to p=pu/2 on Matlock p-y curve
%yt=deformation boundary condition: deflection applied at pile top(m)
%Mt=force boundary condition: bending moment applied at pile top(kN-m)
%ye=elastic solution of y(m), which is used as an initial guess in Matlock's iteration
xpt=zeros(n+5,1);%depth below pile top(m)
x=zeros(n+5,1);%depth below pre-scour ground(m)
for i=3:n+3
    xpt(i)=h*(i-3);
    x(i)=xpt(i)-Lco;
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
itrmax=1000;%maximum number of iterations
tol=1.0e-3;%tolerance for convergence control
y=ye;%use elastic response as the initial guess for y
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
%=================================IN THE END===============================
