
taumax=1; %%% defines T=1
xmin=0; %%% defines the range for space variable
xmax=1;
Nx=15; %%% number of divisions in space
%Ntau=10;
Ntau=6*Nx^2; %%% number of divisions in time
hx=(xmax-xmin)/Nx; %%%% discretization parameter in space
htau=taumax/Ntau; %%%% discretization parameter in time
lambda=htau/hx^2;




solexacta=@(x,t)(exp(-pi^2*t).*sin(pi*x)); %%%% defines the exact solution

ptsx=repmat(linspace(xmin,xmax,Nx+1)',1,Ntau+1); %%%% defines the grid
ptsy=repmat(linspace(0,taumax,Ntau+1),Nx+1,1);
tabsol=solexacta(ptsx,ptsy); %%%% defines the grid function with the exact values

u0=@(x)sin(pi*x); %%%% defines the function for the initial condition
uxmin=@(t)(zeros(size(t))); %%%% defines the functions for the boundary condition
uxmax=@(t)(zeros(size(t)));

solucao=zeros(Nx+1,Ntau+1); %%%% initializes the grid function a zero grid function
solucao(1:Nx+1,1)=u0(linspace(xmin,xmax,Nx+1)); %%%% defines the initial condition
tempos=linspace(0,taumax,Ntau+1);
solucao(1,1:Ntau+1)=uxmin(tempos); %%%% defines the boundary conditions
solucao(end,1:Ntau+1)=uxmax(tempos);

for j=1:Ntau %%%% loop for calculating the solution, by "moving the time forward"
indices=2:Nx;
solucao(indices,j+1)=lambda*(solucao(indices-1,j)+solucao(indices+1,j))+(1-2*lambda)*solucao(indices,j);  
end



figure(1) %%%% plot of the numerical solution
surf(ptsx,ptsy,solucao,'FaceLighting','phong','EdgeColor','none','FaceColor','interp')
title('Approximate solution')
xlabel('x')
ylabel('t')

figure(2) %%%% plot of the absolute error
surf(ptsx,ptsy,abs(tabsol-solucao),'FaceLighting','phong','EdgeColor','none','FaceColor','interp')
title('Absolute error')
xlabel('x')
ylabel('t')



normaerro=max(max(abs(tabsol-solucao))) %%%% norm of the error in infinity-norm














