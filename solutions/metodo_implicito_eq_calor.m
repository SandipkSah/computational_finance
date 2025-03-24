
taumax=1;
xmin=0;
xmax=1;
Ntau=800;
Nx=30;
hx=(xmax-xmin)/Nx;
htau=taumax/Ntau;
lambda=htau/hx^2

solexacta=@(x,t)(exp(-pi^2*t).*sin(pi*x));

ptsx=repmat(linspace(xmin,xmax,Nx+1)',1,Ntau+1);
ptsy=repmat(linspace(0,taumax,Ntau+1),Nx+1,1);
tabsol=solexacta(ptsx,ptsy);

u0=@(x)sin(pi*x);
uxmin=@(t)(zeros(size(t)));
uxmax=@(t)(zeros(size(t)));

solucao=zeros(Nx+1,Ntau+1);
solucao(1:Nx+1,1)=u0(linspace(xmin,xmax,Nx+1));
tempos=linspace(0,taumax,Ntau+1);
solucao(1,1:Ntau+1)=uxmin(tempos);
solucao(end,1:Ntau+1)=uxmax(tempos);



matriz=diag((1+2*lambda)*ones(1,Nx-1))+diag(-lambda*ones(1,Nx-2),1)+diag(-lambda*ones(1,Nx-2),-1);
matriz=sparse(matriz);


for j=1:Ntau
indices=2:Nx;    
vecb=solucao(indices,j);    
vecb(1)=vecb(1)+lambda*solucao(1,j+1);  
vecb(end)=vecb(end)+lambda*solucao(Nx+1,j+1);    
vecx=matriz\vecb;
solucao(indices,j+1)=vecx(:);
end



figure(1)
surf(ptsx,ptsy,solucao,'FaceLighting','phong','EdgeColor','none','FaceColor','interp')
title('Solucao Aproximada')
xlabel('x')
ylabel('t')





figure(2)
surf(ptsx,ptsy,abs(tabsol-solucao),'FaceLighting','phong','EdgeColor','none','FaceColor','interp')
title('Erro absoluto')
xlabel('x')
ylabel('t')





max(max(abs(tabsol-solucao)))



