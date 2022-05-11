figure

% Visualization parameters
N=30; %show last N dimers
sp=round(mean(last_elem))-N; %Visualization starts with (sp) dimer
ep=0; %Visualization ends with (last_elem-ep) dimer
ang=pi/1.2; %rad, projection angle

%Arrays initialization
RxyT=[];
Rxy=[];
Rz=[];
Rz1=[]; Rxy1T=[]; Rxy1D=[]; Rxy_cirD=[]; Rxy_cirT=[];
r=4.8; % nm, Monomer radius
p=3; % helical pitch in monomers
R=19; % nm, MT radius
teta_def=0.2; % rad, angle of monomer bending
phi=linspace(0,2*pi,N_PFs+1);

%calculation of tubulin coordinates
for k = 1:N_PFs
    teta=teta_def;
    tetaz=teta_def;
    nc=0;
    mz=1;
    Rz(k,mz)=0;
    Rxy(k,mz)=R;
    for m = [sp:last_elem(k)-ep]
        for mono=1:2 %monomer cycle
        mz=mz+1;
        nc=nc+1*(C(k,m)==2); % counter for curved dimers
        Rxy(k,mz)=Rxy(k,mz-1)+2*r*sin(teta)*(C(k,m)==2);
        RxyT(k,mz)=Rxy(k,mz)*(G(k,m)>=2);
        Rz(k,mz)=(2*r*mz+2*r*p/13*k)*(C(k,m)==1)+(Rz(k,mz-1)+2*r*cos(tetaz))*(C(k,m)==2);
        teta=teta+0.2*(C(k,m)==2);
        tetaz=tetaz+(0.2-0.05*(mod(nc,5)==0))*(C(k,m)==2); % addition to the angle so that too long ram horns won't overlap
        end
    end
end
        
Rz1=Rz./(Rz(1,3)-Rz(1,2));
Rz1=(Rz1-Rz1(1,2));
Rz1(:,1)=[]; 
Rxy1T=RxyT./19*2.6;
Rxy1D=Rxy./19*2.6;
Rxy1T(:,1)=[];
Rxy1D(:,1)=[];
Rxy1T(Rxy1T==0)=NaN;
Rxy1D(Rxy1D==0)=NaN;

%PFs projections
for n=1:N_PFs
    Rxy_cirT(n,:)=Rxy1T(n,:).*(sin(phi(n)+ang));
    Rxy_cirD(n,:)=Rxy1D(n,:).*(sin(phi(n)+ang));
end

Rz1=4*Rz1; % in nm
Rxy_cirT=4*Rxy_cirT;
Rxy_cirD=4*Rxy_cirD;

% Drawing MT
hold on
for n=[6:12,1:5,13]
    if sum(n==[1:5,13])
        plot(Rxy_cirD(n,:),Rz1(n,:),'-','Color',[0.53,0.86,0.58],'LineWidth',1.5)
        plot(Rxy_cirT(n,:),Rz1(n,:),'-','Color',[1.00,0.45,0.35],'LineWidth',1.5)
    else
        plot(Rxy_cirD(n,:),Rz1(n,:),'-','Color',[0.40,0.6,0.47],'LineWidth',1.5)
        plot(Rxy_cirT(n,:),Rz1(n,:),'-','Color',[0.69,0.26,0.10],'LineWidth',1.5)
    end
end
% Black dots for PF tips
for n=1:13
    indxy=find(~isnan(Rxy_cirD(n,:)),1,'last');
    scatter(Rxy_cirD(n,indxy),Rz1(n,indxy),10,'MarkerEdgeColor',[0.10,0.1,0.10],...
        'MarkerFaceColor',[0.30,0.30,0.30],'MarkerFaceAlpha',0.8,'MarkerEdgeAlpha',0.9)
end

axis image
set(gcf,'Color', [1 1 1])
xlim(4*[-30,30])
ylabel('Length of the visible MT part, nm')