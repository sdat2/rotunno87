%
% This script loads output from the Rotunno-Emanuel (1987) model as encoded
% in hurricane.f, calculates certain additional quantities from that
% output, and graphs the results. 
%
% The user is first prompted for the output file, which should have the
% name 'output_XXX', where XXX is a descriptor (can be any number of
% characters). When prompted, enter only the XXX
%
% The user in the prompted for the number of days into the integration for
% which contour plots are made. (The result is rounded to the nearest
% output time.) Enter 'return' if the last output time is desired. 
%
% Finally, the user is asked for the maximum radius to use in the plots.
% This is helpful to focus on the inner core, for example. Enter an
% arbitrarily large number to plot the whole domain.
%
% Last modified June, 2020
%--------------------------------------------------------------------------
cmap='pink';  % MATLAB colormap for contour plots
%--------------------------------------------------------------------------
%
%  Specify output file to use, output time for contour plots, and outer
%  radius for all plots for which radius is one dimension
%
%sfile=input('Enter output directory (hit return to use "output" )','s');
%if isempty(sfile)
%   sfile='output';
%else
%    sfile=strcat('output_',sfile);
%end  
sfile="output_init";
%
plottime=input('Specifiy time (days) for contour plots (blank for final time)');
%
router=input('Choose outer limit of radial plots (km)');
%
% Open input file 'hurr.in' used to control integration and extract certain
% parameters
%
slash='/';
fid=fopen(strcat(sfile,slash,'hurr.in'),'r');
% Account for addition to hurr.in info on July 1 2020
nlines=16;
fileinfo=dir(strcat(sfile,slash,'hurr.in'));
fileinfo.date % checking value

if datenum(fileinfo.date) > 7.379706613078705e+05
    nlines=17;
end
if datenum(fileinfo.date) > 7.379773968740000e+05    
    nlines=19;
end    
%
C=textscan(fid,'%s %s %*[^\n]',nlines,'headerlines',2);
D=textscan(fid,'%s %s %*[^\n]',5,'headerlines',3);
E=textscan(fid,'%s %s %*[^\n]',14,'headerlines',3);
g=fclose(fid);
%
fields=C{1,1};
pds=str2double(fields(1));
f=1.0e-5*str2double(fields(2));
CD=1.0e-3.*str2double(fields(3));
Ck=1.0e-3.*str2double(fields(4));
pdep=char(fields(9));
diss=str2double(fields(12));
%
fields3=E{1,1};
rwall=str2double(fields3(5));
router=min(router,rwall);
pltime=fields3(10);
endtime=str2double(fields3(8));
%
% Calculate contour file number to use depending on time 
%
clear us usq
dint=str2double(pltime);
ulast=floor(endtime/dint)-1;
days=dint:dint:endtime+dint;
if isempty(plottime) == 0
    [~,idays]=min(abs(plottime-days));
    if idays > ulast
        us='';
        usq=num2str(ulast,'%02d');
    else    
        us=num2str(idays,'%02d');
        usq=us;
    end    
else
    us='';
    usq=num2str(ulast,'%02d');
end    
%
% Open model output files
%
file=strcat(sfile,slash,'time.out');
timetemp=load(file);
time=timetemp;
file=strcat(sfile,slash,'radius.out');
radius=load(file);
rb=radius(:,1);
time1=time(:,1);
file=strcat(sfile,slash,'rgraph.out');
r=load(file);
file=strcat(sfile,slash,'zgraph.out');
z=load(file);
file=strcat(sfile,[slash,'ucon',us,'.out']);
ucon=load(file);
file=strcat(sfile,[slash,'wcon',us,'.out']);
wcon=load(file);
file=strcat(sfile,[slash,'tfcon',us,'.out']);
tfcon=load(file);
%
[nj, mi]=size(tfcon); % Size of most output 2-D files
%
file=strcat(sfile,[slash,'pcon',us,'.out']);
dum=load(file);
%
% Here we account for the fact that the final 2-D pressure field is output
% as a perturbation, whereas otherwise it is output as the full pressure
%
if isempty(us)
    pcon=dum;
else    
    pfull=dum(1:mi,1:nj)'; 
    pconm=mean(pfull,2);
    pconm=repmat(pconm,[1,mi]);
    pcon=pfull-pconm;
end
%
file=strcat(sfile,[slash,'vcon',us,'.out']);
vcon=load(file);
file=strcat(sfile,[slash,'qcon',usq,'.out']);
dum=load(file);
qcon=dum(1:mi,1:nj)';
file=strcat(sfile,[slash,'tescon',us,'.out']);
tescon=load(file);
file=strcat(sfile,[slash,'tecon',us,'.out']);
tecon=load(file);
file=strcat(sfile,[slash,'liqcon',us,'.out']);
liqcon=load(file);
%
% Load input sounding file and extract sea surface temperature from it
%
fid2=fopen(strcat(sfile,slash,'s.in'));
B=textscan(fid2,'%s %s',1);
g2=fclose(fid2);
sst=str2double(B{1,2});
sstk=sst+273.15;
es0=6.112.*exp(17.67.*sst./(243.5+sst)); % Saturation vapor pressure of sea surface
%
dz=z(3)-z(2);
dr=r(3)-r(2);
%
tvcon=tfcon.*(1+0.61*1.0e-3*qcon-0.001*liqcon); % Virtual temperature
% 
% Calculate base state pressure
%
p0=zeros(1,nj);
p0(1)=1000;
for i=2:nj
    p0(i)=p0(i-1).*exp(-9.8.*1e3.*dz./(287.*0.5.*(tvcon(i,mi)+tvcon(i-1,mi))));
end
%
% Calculate full pressure if not using the final output pressure file
%
if isempty(us)
    pfull=zeros(nj,mi);
    pfull(1,1:mi)=p0(1)+pcon(1,1:mi);
    %
    for i=2:nj
        for j=1:mi
            pfull(i,j)=pcon(i,j)+p0(i);
        end    
    end
end    
%
% Calculate a pseudo-density used in mass continuity...product of density
% and virtual potential temperature, from R & E 1987
%
dens=zeros(nj,mi);
dens(1,1:mi)=(100.*pfull(1,1:mi)).^(1-0.287)/287;
for i=2:nj
    for j=1:mi
        dens(i,j)=(100.*pfull(i,j)).^(1-0.287)/287;
    end    
end
%
% Calculate mass streamfunction psi
%
psi=zeros(nj,mi);
for i=1:mi
   for j=2:nj
      psi(j,i)=psi(j-1,i)-1e6.*dz.*(r(i)-0.5*dr).*0.5.*(dens(j,i)+dens(j-1,i)).*...
          ucon(j,i);
   end
   for j=2:nj
       psi(j,i)=psi(j,i)-psi(nj,i).*z(j)./z(nj);  % Force psi to vanish at model top
   end    
end   
%
% Now correct the variable dens to the true density
%
dens=dens.*(100*pfull).^0.287./tvcon;
%
% Use a subroutine at bottom of this file to calculate the "buoyancy
% potential temperature", theta_B, defined as the temperature air would
% have to have at a specified reference pressure (usually 960 hPa)so that
% when lifted reversibly it would have the same density temperature as the
% air at the given level.
%
theta_B=zeros(nj,mi);
for j=1:mi
    ptemp=pfull(:,j)';Ttemp=tfcon(:,j)'-273.15;qtemp=qcon(:,j)';
    temp = thetaB2(ptemp,Ttemp,qtemp);
    theta_B(:,j)=temp;
end
clear temp dum
theta_B=min(theta_B,max(max(theta_B(1:10,:)))); % So that plots are not overwhelmed by stratosphere
%
%  Crudely estimate a potential intensity (not currently used in any
%  graphs)
%
if strcmp(pdep,'y')
    k0=1005*sstk+2.5e6*0.622.*es0./(pfull(1,:)-es0);
else
    k0=zeros(1,mi)+1005*sstk+2.5e6*0.622.*es0./(pds-es0);
end    
eff=100./tfcon(1,:);
if diss > 0.99
    eff=eff.*tfcon(1,:)./(tfcon(1,:)-100);
end    
katm=1005.*tfcon(1,:)+2.5e6.*1e-3*qcon(1,:);
vpotd=sqrt((Ck/CD).*eff.*(k0-katm));

%  Centered differenced version of gradient wind
%
vgradc=zeros(nj,mi);
%
ri2=r+0.5*dr;
b=1000*f*ri2(2:mi-1);
for j=1:nj
    c=-ri2(2:mi-1).*((1./(2*dr*dens(j,2:mi-1))).*100.*(pfull(j,3:mi)-pfull(j,1:mi-2)));
    vgradc(j,2:mi-1)=0.5*(-b+sqrt(max(b.^2-4*c,0)));
end
%
%  Angular momentum
%
r2d=1000*repmat(r,[nj,1]);
M2d=r2d.*(vcon+0.5*f*r2d);
%
%  Plot output via menus
%
timed=time1/24;
colormap(cmap)
menu1=2;
%
while menu1 > 1
    menu1=menu('Choose plot type:','Stop','Evolution','Radial Plots','Contour Plots','Time-radius sections');
    if menu1 == 1
        break
    elseif menu1 == 2
        menu2=2;
        while menu2 > 1
            menu2=menu('Choose plot','Stop','V_max','r_max','u_min','w_max','p_min');
            if menu2 == 1
                break
            elseif menu2 == 2
                v=time(:,2);
                vth=time(:,7);
                vs=time(:,8);
                h=plot(timed,v,timed,vs,timed,vth);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Maximum surface wind speed (m/s)','fontweight','bold')
                legend('Maximum wind','Maximum surface wind','Potential intensity','location','best') 
                title('Azimuthal Velocities (m/s)')
            elseif menu2 == 3
                rmax=time(:,6);
                h=plot(timed,rmax);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Radius of maximum winds (km)','fontweight','bold')
                title('Radius of Maximum Winds (km)')
            elseif menu2 == 3
                ut=time(:,3);
                h=plot(timed,ut);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Minimum surface radial wind speed (m/s)','fontweight','bold')
                title('Surface Radial Velocity (m/s)')
            elseif menu2 == 4    
                ut=time(:,3);
                h=plot(timed,ut);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Minimum surface radial wind speed (m/s)','fontweight','bold')
                title('Minimum Radial Velocity in Outflow (m/s)')
            elseif menu2 == 5    
                wt=time(:,4);
                h=plot(timed,wt);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Maximum vertical velocity (cm/s)','fontweight','bold')
                title('Maximum Vertical Velocity (cm/s)')
            elseif menu2 == 6
                pc=time(:,5);
                h=plot(timed,pc);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Time (days)','fontweight','bold')         
                ylabel('Central surface pressure (mb)','fontweight','bold')
                title('Central Surface Pressure (mb)')
            end    
        end
        %
    elseif menu1 == 3
        menu3=2;
        while menu3 > 1
            menu3=menu('Choose Plot','Stop','Surface azimuthal velocity', ...
                'Surface radial velocity','Surface pressure','Surface theta_e and mid-level theta_e*');
            if menu3 == 1
                break
            elseif menu3 == 2
                lev=1;
                v=vcon(lev,:);
                rbplus=rb+0.5*dr;
                M=rbplus.*v+0.5*f*rbplus.^2;
                h=plot(ri2,v,'b',ri2,vgradc(lev,:),'g');
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel('Surface azimuthal velocity (m/s)','fontweight','bold')
                legend('Azimuthal wind','Gradient wind')
                title('Surface Azimuthal and Gradient Winds (m/s)')
            elseif menu3 == 3
                ub=ucon(1,:);
                h=plot(rb,ub);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel('Surface radial velocity (m/s)','fontweight','bold')
                title('Surface Radial Velocity (m/s)')
            elseif menu3 == 4
                p=pcon(1,:);
                h=plot(rb,p);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel('Surface pressure (mb)','fontweight','bold')
                title('Surface Pressure (mb)')
            elseif menu3 == 5
                lev1=2;lev2=11;
                file2=strcat(sfile,'/tescon.out');
                tescon=load(file2);
                thes(1,:)=tescon(lev2,:);
                the=tecon(lev1,:);
                h=plot(rb,the,rb,thes);
                set(h,'linewidth',2)
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel('Surface theta_e (K) and mid-level theta_e* (K)','fontweight','bold')
                entry1=['\theta_e at ', num2str(z(lev1)),' km'];
                entry2=['\theta_{es} at ', num2str(z(lev2)),' km'];
                legend(entry1,entry2,'location','best')
                title('Surface theta_e (K) and Mid-level theta_e* (K)')
            end
            set(gca,'xlim',[0 router])
        end
        %
    elseif menu1 == 4
        menu4=2;
        while menu4 > 1
            menu4=menu('Choose Plot','Stop','Azimuthal velocity', ...
                'Radial velocity','Vertical velocity','Streamfunction', ...
                'Perturbation pressure','Perturbation temperature', ...
                'Theta_B','Theta_e','Theta_es','Theta_es and M', ...
                'Liquid water', 'Vertical turbulent diffusivity (final time only)');
            if menu4 ==1
                break
            elseif menu4 == 2
                rm=1000.*r;
                q=max(size(z));
                contourf(r,z,vcon,30)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Azimuthal Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 3
                contourf(r,z,ucon,15)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Radial Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 4
                % Make more contours for negative values
                tempcon=zeros(nj,mi)+1;
                tempcon(wcon<0)=10;
                wconp=wcon.*tempcon;
                contourf(r,z,wconp)
                g=colorbar;
                % Make colorbar labels consistent with contours
                u=get(g,'TickLabels');
                un=str2double(u);
                qun=max(size(un));
                un=un';
                fac=zeros(1,qun)+1;
                fac(un<0)=0.1;
                un=un.*fac;
                set(g,'TickLabels',num2str(un'))
                set(g,'fontweight','bold')
                title('Vertical Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 5
                [~,imax]=min(abs(r-router));
                imax=imax+1;
                contourf(r(1:imax),z,psi(:,1:imax),15)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Mass Streamfunction (Kg/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')    
            elseif menu4 == 6
                [mtot,ntot]=size(pcon);
                pconm=mean(pcon,1);
                pconm=repmat(pconm,[mtot,1]);
                pconp=pcon-pconm;
                pcolor(r,z,pconp(1:nj,1:mi))
                shading interp
                hold on
                contour(r,z,pconp(1:nj,1:mi),'k')
                hold off
                g=colorbar;
                set(g,'fontweight','bold')
                title('Perturbation Pressure (mb)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 7    
                file=strcat(sfile,['/tcon',us,'.out']);
                tcon=load(file);
                if isempty(us) == 0
                    tcon=tcon(1:mi,1:nj)';
                    tconm=repmat(mean(tcon,2),[1,mi]);
                    tcon=tcon-tconm;
                end    
                contourf(r,z,tcon(1:nj,1:mi),10)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Perturbation Temperature (K)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 8     
                contourf(r,z,theta_B,20)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Buoyancy Potential Temperature (K)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')    
            elseif menu4 == 9     
                temax=max(tecon(1,:))+5;
                teconp=min(tecon,temax);
                contourf(r,z,teconp,20)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Equivalent Potential Temperature (K)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')  
            elseif menu4 == 10     
                file=strcat(sfile,'/tecon.out');
                tecon=load(file);
                temax=max(tecon(1,:));
                teconp=min(tecon,temax);
                tesconp=min(tescon,max(teconp(1,:))+3);
                contourf(r,z,tesconp,20)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Saturation Equivalent Potential Temperature (K)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')   
            elseif menu4 == 11     
                %temax=max(tecon(1,:));
                %tecon=min(tecon,temax);
                %contourf(r,z,theta_B,20)
                %caxis([min(min(theta_B)) max(max(theta_B))]);
                tesconp=min(tescon,max(tecon(1,:))+3);
                testemp=tesconp;
                testemp(1:8,:)=tecon(1:8,:); % Replace theta_es by theta_e in boundary layer
                contourf(r,z,testemp,50)
                caxis([min(min(testemp)) max(max(testemp))]);
                %
                hold on
                contour(r,z,sqrt(max(M2d,0)),100,'color','k','linewidth',1.0)
                hold off
                g=colorbar;
                set(g,'fontweight','bold')
                title('\theta_{es} (shading) and Angular Momentum (black)','fontweight','bold')
                set(gca,'fontweight','bold','ylim',[0,17])
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')     
            elseif menu4 == 12     
                %liqcon=log(liqcon+0.1);   
                contourf(r,z,liqcon,20)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Liquid Water Content (g Kg^{-1})','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')
            elseif menu4 == 13     
                file=strcat(sfile,'/xkcon.out');
                xkcon=load(file);
                contourf(r,z,xkcon,10)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Vertical turbulent diffusivity (m^2/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('z (km)','fontweight','bold')               
            end
            set(gca,'xlim',[0 router])
        end
        %
    elseif menu1 == 5
        menu5=2;
        while menu5 > 1
            menu5=menu('Choose Plot','Stop','Surface azimuthal velocity', ...
                'Surface radial velocity','Surface theta_e', ...
                '15 km azimuthal velocity', '15 km radial velocity');
            if menu5 == 1
                break
            elseif menu5 == 2
                file=strcat(sfile,'/vbhov.out');
                vhov=load(file);
                contourf(r,timed,vhov,10)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Azimuthal Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('Time (days)','fontweight','bold')
            elseif menu5 == 3
                file=strcat(sfile,'/ubhov.out');
                uhov=load(file);
                contourf(r,timed,uhov,10)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Radial Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('Time (days)','fontweight','bold')    
            elseif menu5 == 4    
                file=strcat(sfile,'/thehov.out');
                thehov=load(file);
                contourf(r,timed,thehov,10)
                g=colorbar;
                set(g,'fontweight','bold')
                title('Surface theta_e (K)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('Time (days)','fontweight','bold')
            elseif menu5 == 5
                file=strcat(sfile,'/vthov.out');
                vthov=load(file);
                contourf(r,timed,vthov,10)
                title('15 km Azimuthal Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('Time (days)','fontweight','bold')
                g=colorbar;
                set(g,'fontweight','bold')
            elseif menu5 == 6
                file=strcat(sfile,'/uthov.out');
                uthov=load(file);
                contourf(r,timed,uthov,10)
                title('15 km Radial Velocity (m/s)','fontweight','bold')
                set(gca,'fontweight','bold')
                xlabel('Radius (km)','fontweight','bold')
                ylabel ('Time (days)','fontweight','bold')
                g=colorbar;
                set(g,'fontweight','bold')                
            end
            set(gca,'xlim',[0 router])
        end
    end
end    
%
function [ theta_B ] = thetaB2(p,T,q )
% This function takes vectors of pressure (hPa), temperature (C), and
% mixing ratio (g/Kg) and calculates the buoyancy potential temperature,
% theta_B (K).
%
% October, 2016. Modified June, 2020 to use partial pressure of dry air
%--------------------------------------------------------------------------
% Constants
%
pnorm=960; % Pressure level at which theta_B defined
tolerance=0.02; % Accuracy of theta_B (K)
%
cp=1005;
cpv=1870;
Rd=287;
Cl=4160;  
Rv=461;
Lv0=2.5e6;
pref=1000;
Tref=290;
%
Tk=T+273.15;
qq=0.001.*q;
n=max(size(p));
theta_B=zeros(1,n);
%
Tv=Tk.*(1+qq./0.622)./(1+qq)-273.15; % Virtual temperature
Lv=Lv0+(cpv-Cl).*T;
es=6.112*exp(16.67.*T./(243.5+T));
qs=(Rd/Rv).*es./(p-es);
s=(cp+Cl.*qs).*log(Tk./Tref)-Rd.*log((p-es)./pref)+Lv.*qs./Tk; % Saturation entropy
%
[Tg,qsg]=ssinvert(s,pnorm); % Initial guess of theta_B
%
ibot=1e6;
for i=1:n
    if p(i) < pnorm
        ibot=min(i,ibot);
        del=10;
        while abs(del) > tolerance
            Lv=Lv0+(cpv-Cl).*(Tg(i)-273.15);
            sg=(cp+Cl*qsg(i))*log(Tg(i)/Tref)-Rd*log((pnorm-es(i))/pref)+Lv*qsg(i)./Tg(i);
            [Tnew,qsnew]=ssinvert(sg,p(i));
            Tvnew=Tnew.*(1+qsnew./0.622)./(1+qsg(i))-273.15;
            del=Tv(i)-Tvnew;
            Tg(i)=Tg(i)+0.2*del;
            Tgc=Tg(i)-273.15;
            es(i)=6.112*exp(16.67.*Tgc./(243.5+Tgc));
            qsg(i)=(Rd/Rv)*es(i)./(pnorm-es(i));
        end
    %    
    end    
    theta_B(i)=Tg(i);
end
%
for i=1:ibot
    theta_B(i)=theta_B(ibot);
end    
%}
end
%
function [ T,qs] = ssinvert( ss,p )
%  This function accepts a vector of parcel saturation entropy, ss,
%  and a scalar pressure p (hPa), and returns vectors of temperature (K),
%  and saturation specific humidity (g/g) along a reversible adiabat.
%  Newton-Raphson method.
%  October, 2016. Modified June, 2020 to use partial pressure of dry air
%
cp=1005;
cpv=1870;
Rd=287;
Cl=4160;  
Rv=461;
Lv0=2.5e6;
pref=1000;
Tref=290;
tolerance=0.02;
%
sssize=max(size(ss));
T=zeros(1,sssize)+280;
qs=zeros(1,sssize);
%
for i=1:sssize
    del=1;
    while abs(del) > tolerance 
        %
        Tc=T(i)-273.15;
        Lv=Lv0+(cpv-Cl)*Tc;
        esg=6.112*exp(17.67*Tc/(243.5+Tc));
        qs(i)=0.622*esg/(p-esg);
        snew=(cp+Cl*qs(i))*log(T(i)/Tref)-Rd*log((p-esg)/pref)+Lv*qs(i)/T(i);
        dsdT=(cp+Cl*qs(i))/T(i)+Lv^2*qs(i)/(Rv*T(i)^3);
        dsdT=dsdT+Rd*Lv*esg/(Rv*T(i)^2*(p-esg));
        del=ss(i)-snew;
        T(i)=T(i)+del/dsdT;
    end
    Tc=T(i)-273.15;
    esg=6.112*exp(17.67*Tc/(243.5+Tc));
    qs(i)=0.622*esg/(p-esg);
end
end