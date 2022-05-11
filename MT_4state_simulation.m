clear

% Define microtubule (MT) dimensions and preallocate memory for simulation
N_PFs=13; % Number of protofilaments (PFs) in microtubule
N_layers=10000; % preallocated space for maximal MT length in dimers (=80 um)
C=zeros(N_PFs,N_layers); % array for tubulin conformations (0-no tubulin, 1-straight, 2-curved)
G=C; % array for tubulin nucleotide states (0-no tubulin, 1-gdp, 2-gtp)
GC=C; % tubulin nucleotide states without curved dimers (only those tubulins that contribute to lateral energy)
LatBonds=C; % lateral energy for bending calculation
LongBonds=C; % longitudinal energy for detachment calculation
dGlat=zeros(4,4); % lateral bonds energies for different nucletide combinations

% Technical parameters
calc_depth=4; % depth of lateral energy matrix update
npass=800; % Output variable (time, mtlength, etc) will store every npass frame of the simulation
% (high npass values speed up calculations and free up RAM, but cause low temporal resolution in the output trajectory) 
NUMBER_OF_TIME_STEPS=500*900; % number of planned iterations (assuming npass is 800, 400 sparsed iterations = 1 min of simulation)

% Define initial MT configuration
Seed_length=20; % Length of a superstable MT seed in dimers (will not disassemble after catastrophe)
G(:,[1:Seed_length])=3*ones(N_PFs,Seed_length); % seed nucleotide (3 - is a state of a special superstable nucleotide)
C(:,[1:Seed_length])=ones(N_PFs,Seed_length); % seed conformation (1 - is a straight state)
GC(C==1)=G(C==1); % only straight dimers

% Define simulation parameters
calc_time=60; % s, Time to calculate dynamics
ctub=10; % tubulin concentration, uM
dGlong(1)=-12.8; dGlong(2)=-13.2; %kT, longitudinal energy bond (1 - GDP, 2 - GTP)
kplus=0.55; %1/(uM*s), association constant
Ebend(1)=8; Ebend(2)=6; Ebend(3)=6.5; % kT, bending energy (1 - GDP, 2 - GTP, 3 - seed)
dGlat(2,2)=-5.6; %kT, lateral energy bond between GDP dimers
dGlat(3,3)=-7.5; %kT lateral energy bond between GTP dimers
dGlat([7;10])=-7.5; %kT lateral energy bond between GTP and GDP dimers ([3,2] & [2,3])
dGlat(4,2:4)=-20; dGlat(2:4,4)=-20; % kT, lateral energy bond between dimers in the seed
khydr=0.09; %1/s hydrolysis constant
kstr_gtp=300; %1/s, straightening contant for GTP-tubulin interface
kstr_gdp=50; %1/s, straightening contant for GDP-tubulin interface
lambda=105; %barrier coefficient
% aggregating all kstr in one array for convenience
kstr(1,1)=kstr_gdp; kstr(1,2)=kstr_gdp; kstr(1,3)=kstr_gdp/lambda; %s-1, GDP tubulins (1-zero lateral bonds, 2-one lateral bond, 3-two lateral bonds)
kstr(2,1)=kstr_gtp; kstr(2,2)=kstr_gtp; kstr(2,3)=kstr_gtp/lambda; %s-1, GTP tubulins -//-
kstr(3,1)=100; kstr(3,2)=100; kstr(3,3)=100; %s-1, Seed -//- 

% Lateral bonds calculation
Latf=@(a,b) dGlat(((b+1)-1)*size(dGlat,1)+(a+1)); % function to calculate lateral energy (linear indexing)
LatBonds(1,2:end)=Latf(GC(1,2:end),GC(2,2:end)); % at the seam, right
LatBonds(1,2:end)=LatBonds(1,2:end)+Latf(GC(1,2:end),GC(N_PFs,1:end-1))/2; % at the seam, left top
LatBonds(1,3:end)=LatBonds(1,3:end)+Latf(GC(1,3:end),GC(N_PFs,1:end-2))/2; % at the seam, left bottom
LatBonds(2:N_PFs-1,:)=reshape(Latf(GC(2:N_PFs-1,:),GC(1:N_PFs-2,:)),N_PFs-2,N_layers); % to the left of the seam
LatBonds(2:N_PFs-1,:)=LatBonds(2:N_PFs-1,:)+reshape(Latf(GC(2:N_PFs-1,:),GC(3:N_PFs,:)),N_PFs-2,N_layers); % to the right of the seam
LatBonds(N_PFs,:)=Latf(GC(N_PFs,:),GC(N_PFs-1,:)); % at the seam, left
LatBonds(N_PFs,1:end-1)=LatBonds(N_PFs,1:end-1)+Latf(GC(N_PFs,1:end-1),GC(1,2:end))/2; % at the seam, right bottom
LatBonds(N_PFs,1:end-2)=LatBonds(N_PFs,1:end-2)+Latf(GC(N_PFs,1:end-2),GC(1,3:end))/2; % at the seam, right top
% Longitudinal bonds calculations
LongBonds(:,2:end)=1*(G(:,1:end-1)==1)+2*(G(:,1:end-1)>1);

% initialize arrays of MT characteristics
time=zeros(1,NUMBER_OF_TIME_STEPS);
mtlength=time; %Medium length of straight PFs parts
ramhorns=time; %Mean curved PFs length
sdstraight=time; %Raggedness (standard deviation) of straight PFs parts
gtp_cap=zeros(N_PFs,NUMBER_OF_TIME_STEPS); %Number of GTP-tubulin dimers for each PF
pflength=gtp_cap; %Lengths of straight parts for each PF
Tot_energy=zeros(1,N_PFs); %Total energy of energy needed to bend dimer
last_elem=zeros(1,N_PFs); %last tubulin dimer in a PF
ind_lat=zeros(N_PFs,1); %number of lateral bonds for bending calculation
for le=1:N_PFs
last_elem(le)=find(C(le,:),1,'last'); 
end
last_straight_col(1:N_PFs,1)=last_elem; % last dimer in a straight state of each PF
curv_row=[]; % indices of a PF where straightening may take place
curv_col=[]; % indices of a dimer in a PF -//-
hyd_row=[]; % indices of a PF where GTP hydrolysis may take place
hyd_col=[]; % indices of a dimer in a PF -//-

time(1)=0;
temp_time=0;
mtlength(1)=Seed_length;
it=1;
itcut=it; % sparsed iterations

while time(itcut) < calc_time
    it=it+1;
        
    if isempty(curv_row)
        curv_row=[];
        curv_col=[];
    end
    if isempty(hyd_row)
        hyd_row=[];
        hyd_col=[];
    end
    
    % determine the times of possible events
    
%     attachment
    [tatt,att_j]=min(-log(rand(1,N_PFs))./(kplus*ctub));
    
    tdet=inf;
    tstr=inf;
    if ~isempty(curv_row)
        % detachment
        i_long=LongBonds((curv_col-1)*size(C,1)+curv_row); % determining the nucleotide of longitudinal bond to break
        koff=kplus*10^6./exp(-dGlong(i_long));
        [tdet,i_det]=min(-log(rand(1,length(curv_row)))./koff);
        det_j=curv_row(i_det); % PF of detachment
        det_i=curv_col(i_det); % position of detachment
        
        % straightening
        un=sort(curv_row);
        un(un((1:end-1))==un(2:end))=[]; % Choose only unique PFs
        
        ind_lats=0;
        % determine how many lateral neighbors the PF will have if it is straightened
        for ii = 1:length(un)
            uni=un(ii);
            if uni==1 % at the seam
                ind_lats(ii)=1*(GC(N_PFs,last_straight_col(uni)-1)==0 && GC(uni+1,last_straight_col(uni)+1)==0)+...
                    2*(GC(N_PFs,last_straight_col(uni)-1)~=0 && GC(uni+1,last_straight_col(uni)+1)==0)+...
                    2*(GC(N_PFs,last_straight_col(uni)-1)==0 && GC(uni+1,last_straight_col(uni)+1)~=0)+...
                    3*(GC(N_PFs,last_straight_col(uni)-1)~=0 && GC(uni+1,last_straight_col(uni)+1)~=0);
            elseif uni==N_PFs % at the seam
                ind_lats(ii)=1*(GC(uni-1,last_straight_col(uni)+1)==0 && GC(1,last_straight_col(uni)+3)==0)+...
                    2*(GC(uni-1,last_straight_col(uni)+1)~=0 && GC(1,last_straight_col(uni)+3)==0)+...
                    2*(GC(uni-1,last_straight_col(uni)+1)==0 && GC(1,last_straight_col(uni)+3)~=0)+...
                    3*(GC(uni-1,last_straight_col(uni)+1)~=0 && GC(1,last_straight_col(uni)+3)~=0);
            else % at any other PF
                ind_lats(ii)=1*(GC(uni-1,last_straight_col(uni)+1)==0 && GC(uni+1,last_straight_col(uni)+1)==0)+...
                    2*(GC(uni-1,last_straight_col(uni)+1)~=0 && GC(uni+1,last_straight_col(uni)+1)==0)+...
                    2*(GC(uni-1,last_straight_col(uni)+1)==0 && GC(uni+1,last_straight_col(uni)+1)~=0)+...
                    3*(GC(uni-1,last_straight_col(uni)+1)~=0 && GC(uni+1,last_straight_col(uni)+1)~=0);
            end
        end
        
        lin_nuc=(last_straight_col(un)-1)*N_PFs+un;
        nucl=G(lin_nuc);
        
        [tstr,j]=min(-log(rand(1,length(un)))./kstr((ind_lats-1)*3+nucl'));
        str_j=un(j);
    end
    
    % bending
    lin_nuc=(last_straight_col-2)*N_PFs+(1:N_PFs)';
    nucl=G(lin_nuc); %nucleotide state of the interface to be bent
    Tot_energy = LatBonds((last_straight_col-1)'*N_PFs+(1:N_PFs))+Ebend(nucl); %Total energy needed to bend each PF
    
    % determine how many lateral bonds the PF needs to break
    for ii = 1:N_PFs
        if ii==1 % at the seam
            ind_lat(ii,1)=1*(GC(N_PFs,last_straight_col(ii)-2)==0 && GC(ii+1,last_straight_col(ii))==0)+...
                2*(GC(N_PFs,last_straight_col(ii)-2)~=0 && GC(ii+1,last_straight_col(ii))==0)+...
                2*(GC(N_PFs,last_straight_col(ii)-2)==0 && GC(ii+1,last_straight_col(ii))~=0)+...
                3*(GC(N_PFs,last_straight_col(ii)-2)~=0 && GC(ii+1,last_straight_col(ii))~=0);
        elseif ii==N_PFs % at the seam
            ind_lat(ii,1)=1*(GC(ii-1,last_straight_col(ii))==0 && GC(1,last_straight_col(ii)+2)==0)+...
                2*(GC(ii-1,last_straight_col(ii))~=0 && GC(1,last_straight_col(ii)+2)==0)+...
                2*(GC(ii-1,last_straight_col(ii))==0 && GC(1,last_straight_col(ii)+2)~=0)+...
                3*(GC(ii-1,last_straight_col(ii))~=0 && GC(1,last_straight_col(ii)+2)~=0);
        else % at any other PF
            ind_lat(ii,1)=1*(GC(ii-1,last_straight_col(ii))==0 && GC(ii+1,last_straight_col(ii))==0)+...
                2*(GC(ii-1,last_straight_col(ii))~=0 && GC(ii+1,last_straight_col(ii))==0)+...
                2*(GC(ii-1,last_straight_col(ii))==0 && GC(ii+1,last_straight_col(ii))~=0)+...
                3*(GC(ii-1,last_straight_col(ii))~=0 && GC(ii+1,last_straight_col(ii))~=0);
        end
    end
    
    [tbend,bend_j]=min(-log(rand(1,N_PFs)).*exp(-Tot_energy)./kstr((ind_lat-1)*3+nucl)');
    
    thydr=inf;
    % hydrolysis
    if ~isempty(hyd_row)
        [thydr,i_hydr]=min(-log(rand(1,length(hyd_row)))./khydr);
        hydr_j=hyd_row(i_hydr); % PF of hydrolysis
        hydr_i=hyd_col(i_hydr); % position of hydrolysis
    end
    
    [min_t, event]=min([tatt tdet tstr tbend thydr]); % event selection
    
    if event==1 % attachment
        att_i=last_elem(att_j)+1;
        C(att_j,att_i)=2; % new dimer is in curved state
        G(att_j,att_i)=2; % GTP-state
        last_elem(att_j)=last_elem(att_j)+1;
        curv_row(end+1,1)=att_j;
        curv_col(end+1,1)=att_i;
        ch_d=att_i-calc_depth;
        ch_j=att_j;
    end
    
    if event==2 % detachment
        det_ran=det_i:last_elem(det_j);
        C(det_j,det_ran)=0;
        G(det_j,det_ran)=0;
        GC(det_j,det_ran)=0;
        tempi=curv_row==det_j & sum(curv_col==det_ran,2);
        curv_row(tempi)=[];
        curv_col(tempi)=[];
        if ~isempty(hyd_row) % forget potential hydrolysis sites in detached PF
            temph=hyd_row==det_j & sum(hyd_col==det_i-1:last_elem(det_j),2);
            hyd_row(temph)=[];
            hyd_col(temph)=[];
        end
        last_elem(det_j)=last_elem(det_j)-(last_elem(det_j)-det_i+1);
        ch_d=det_i-calc_depth;
        ch_j=det_j;
    end
    
    if event==3 % straightening
        str_i=min(curv_col(curv_row==str_j));
        C(str_j,str_i)=1; % straight
        GC(str_j,str_i)=G(str_j,str_i);
        last_straight_col(str_j,1)=last_straight_col(str_j,1)+1;
        tempi=curv_row==str_j & curv_col==str_i;
        curv_row(tempi)=[];
        curv_col(tempi)=[];
        if G(str_j,str_i-1)==2 
            if ~isempty(hyd_row) % add a new potential hydrolysis site for a dimer under the straight one
                temph=hyd_row==str_j & hyd_col==str_i-1;
                if sum(temph)==0
                    hyd_row(end+1,1)=str_j;
                    hyd_col(end+1,1)=str_i-1;
                end
            else
                hyd_row(end+1,1)=str_j;
                hyd_col(end+1,1)=str_i-1;
            end
        end
        ch_d=str_i-calc_depth;
        ch_j=str_j;
    end
    
    if event==4 %bending
        tempi=last_straight_col(bend_j);
        C(bend_j,tempi)=2; % bent state
        GC(bend_j,tempi)=0;
        last_straight_col(bend_j,1)=last_straight_col(bend_j)-1;
        curv_row(end+1,1)=bend_j;
        curv_col(end+1,1)=tempi;     
        if ~isempty(hyd_row)
            temph=hyd_row==bend_j & hyd_col==tempi-1;
            hyd_row(temph)=[];
            hyd_col(temph)=[];
        end
        ch_d=last_straight_col(bend_j)-calc_depth;
        ch_j=bend_j;
    end
    
    if event==5 %hydrolysis
        G(hydr_j,hydr_i)=1; %GTP->GDP
        GC(hydr_j,hydr_i)=C(hydr_j,hydr_i)==1;
        hyd_row(i_hydr)=[];
        hyd_col(i_hydr)=[];
        ch_d=hydr_i-calc_depth;
        ch_j=hydr_j;
    end
    
    % renewal of energies in a region near the changes caused by last accepted event
    end_opt=max(last_elem);
    range=ch_d:end_opt;
    if ch_j == 1
        LatBonds(1,range)=Latf(GC(1,range),GC(2,range)); % at the seam, right
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-2))/2; % at the seam, left bottom
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-1))/2; % at the seam, left top
        LatBonds(2,range)=Latf(GC(2,range),GC(1,range)); % to the left of the seam
        LatBonds(2,range)=LatBonds(2,range)+Latf(GC(2,range),GC(3,range)); % to the right of the seam
        LatBonds(N_PFs,range)=Latf(GC(N_PFs,range),GC(N_PFs-1,range)); % at the seam, left
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+1))/2; % at the seam, right bottom
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+2))/2; % at the seam, right top
    elseif ch_j == 2
        LatBonds(1,range)=Latf(GC(1,range),GC(2,range)); % at the seam, right
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-2))/2; % at the seam, left bottom
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-1))/2; % at the seam, left top
        LatBonds(2,range)=Latf(GC(2,range),GC(1,range)); % to the left of the seam
        LatBonds(2,range)=LatBonds(2,range)+Latf(GC(2,range),GC(3,range)); % to the right of the seam
        LatBonds(3,range)=Latf(GC(3,range),GC(2,range)); % to the left of the seam
        LatBonds(3,range)=LatBonds(3,range)+Latf(GC(3,range),GC(4,range)); % to the right of the seam
    elseif ch_j == N_PFs-1
        LatBonds(N_PFs-2,range)=Latf(GC(N_PFs-2,range),GC(N_PFs-3,range)); % to the left of the seam
        LatBonds(N_PFs-2,range)=LatBonds(N_PFs-2,range)+Latf(GC(N_PFs-2,range),GC(N_PFs-1,range)); % to the right of the seam
        LatBonds(N_PFs-1,range)=Latf(GC(N_PFs-1,range),GC(N_PFs-2,range)); % to the left of the seam
        LatBonds(N_PFs-1,range)=LatBonds(N_PFs-1,range)+Latf(GC(N_PFs-1,range),GC(N_PFs,range)); % to the right of the seam
        LatBonds(N_PFs,range)=Latf(GC(N_PFs,range),GC(N_PFs-1,range)); % on seam, left
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+1))/2; % at the seam, right bottom
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+2))/2; % at the seam, right top
    elseif ch_j == N_PFs
        LatBonds(N_PFs-1,range)=Latf(GC(N_PFs-1,range),GC(N_PFs-2,range)); % to the left of the seam
        LatBonds(N_PFs-1,range)=LatBonds(N_PFs-1,range)+Latf(GC(N_PFs-1,range),GC(N_PFs,range)); % to the right of the seam
        LatBonds(N_PFs,range)=Latf(GC(N_PFs,range),GC(N_PFs-1,range)); % at the seam, left
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+1))/2; % at the seam, right bottom
        LatBonds(N_PFs,range)=LatBonds(N_PFs,range)+Latf(GC(N_PFs,range),GC(1,range+2))/2; % at the seam, right top
        LatBonds(1,range)=Latf(GC(1,range),GC(2,range)); % at the seam, right
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-2))/2; % at the seam, left bottom
        LatBonds(1,range)=LatBonds(1,range)+Latf(GC(1,range),GC(N_PFs,range-1))/2; % at the seam, left top
    else
        LatBonds(ch_j-1,range)=Latf(GC(ch_j-1,range),GC(ch_j-2,range)); % to the left of the seam
        LatBonds(ch_j-1,range)=LatBonds(ch_j-1,range)+Latf(GC(ch_j-1,range),GC(ch_j,range)); % to the right of the seam
        LatBonds(ch_j,range)=Latf(GC(ch_j,range),GC(ch_j-1,range)); % to the left of the seam
        LatBonds(ch_j,range)=LatBonds(ch_j,range)+Latf(GC(ch_j,range),GC(ch_j+1,range)); % to the right of the seam
        LatBonds(ch_j+1,range)=Latf(GC(ch_j+1,range),GC(ch_j,range)); % to the left of the seam
        LatBonds(ch_j+1,range)=LatBonds(ch_j+1,range)+Latf(GC(ch_j+1,range),GC(ch_j+2,range)); % to the right of the seam
    end
    
    LongBonds(:,range)=1*(G(:,range-1)==1)+2*(G(:,range-1)>1);
    
    temp_time=temp_time+min_t;
    
    % calculating time and other characteristics of MT every npass steps
    if mod(it,npass)==0 || it==2
        
        itcut=itcut+1;
        ramhorns(itcut)=mean(last_elem'-last_straight_col);
        mtlength(itcut)=median(last_straight_col);
        sdstraight(itcut)=std(last_straight_col);
        
        for nn=1:N_PFs
            gtp_cap(nn,itcut)=sum(G(nn,:)==2);
            pflength(nn,itcut)=last_straight_col(nn);
        end
        
        time(itcut)=time(itcut-1)+temp_time;
        temp_time=0;
    end
    
end

% squeeze arrays
num=1:itcut;
time=time(num);
mtlength=mtlength(num);
ramhorns=ramhorns(num);
sdstraight=sdstraight(num);
gtp_cap=gtp_cap(:,num);
pflength=pflength(:,num);

% outputs
velocity=(mtlength(itcut)-mtlength(1))/time(itcut)*8 % MT growth speed, nm/s
Curved_PFs_length_mean=8*mean(ramhorns) % Curved PF length, nm
Raggednes_mean=8*mean(sdstraight) % Raggedness, nm (standard deviation of straight PF)
GTPcap=sum(sum(G==2)) % GTP content

figure
subplot(3,1,1)
plot(time,8*mtlength,'k-');
xlabel('time, sec');
ylabel('MT length, nm')
subplot(3,1,2)
plot(time,8*ramhorns,'k-');
xlabel('time, sec');
ylabel('Curved PF length,nm')
subplot(3,1,3)
plot(time,8*sdstraight,'k-');
xlabel('time, sec');
ylabel('Raggedness, nm');
figure
plot(time,gtp_cap);
xlabel('time, sec');
ylabel('Total number of GTP-tubulins in each PF');
figure
plot(time,8*pflength)
xlabel('time, sec');
ylabel('PF Length, nm');


