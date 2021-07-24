%load 'MODELO_MU'

dx=25;
[vel,x,z]=marmousimodel(dx);
xmax=max(x);
zmax=max(z); 
vlow=min(min(vel));
vhigh=max(max(vel));
dt=.004;
dtstep=.001; 
tmax=2*zmax/vlow; 
[seisfilt,seis,t]=afd_explode(dx,dtstep,dt,tmax,vel,x,zeros(size(x)),[5 10 40 50],0,2); 
dtau=8*dt;
[vrmsmod,tau]=vzmod2vrmsmod(vel,z,dtau,tmax);
frange=[2 80];
taucheck=0:.1:2;
[zosmigt,exzos]=pspi_stack_tmig_rms(seis,t,x,vrmsmod,x,tau,frange,taucheck); 
zcheck=0:100:2000;
frange=[2 80];
[zosmigd,exzos]=pspi_stack(seis,t,x,vel,x,z,frange,zcheck); 
%%
%[tw, w] = ricker(25, dt);










%%
for i=1:length(x)
    
imp(:, i)=blimp(zosmigd(:,i),vel(:,i), T,15,100);

end


for i=1:length(x)
    Imp(:,i)=seisinv1(zosmigt(:,i),t,vel(1,1),15,100,300);
end 


for i=1:length(x)
[I_recursiva(:,i), s_recursiva(:,i)]=invrecursiva(zosmigd(:,i),vel(:,i),T,20,100, w);
end

for i=1:length(x)
[I_wagn(:,i), s_wagn(:,i)]=invwagn(zosmigd(:,i),T,20,100,w, 0.8, vel(:,i), 5);
end 