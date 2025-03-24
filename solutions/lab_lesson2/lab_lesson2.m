%%%% exercise 1

f=@(x)sin(x);
%x0=0;
x0=0.5; %%%% try these two different choices for x0 and see the differences at the estimates for the order of convergences
exact_value=cos(x0);

hs=1./(10:10:100);

ffd=forward(f,x0,hs); %%% approximations obtained using forward finite difference
bfd=backward(f,x0,hs); %%% approximations obtained using backward finite difference
cfd=centered(f,x0,hs); %%% approximations obtained using centered finite difference

abs_error_f=abs(ffd-exact_value); %%% absolute error of forward finite difference
abs_error_b=abs(bfd-exact_value); %%% absolute error of backward finite difference
abs_error_c=abs(cfd-exact_value); %%% absolute error of centered finite difference

%%%% estimate of the order of convergence - forward f.d.
p=polyfit(log(hs),log(abs_error_f),1);
p(1) %%%% the slope is the estimate for the order of convergence

%%%% estimate of the order of convergence - backward f.d.
p=polyfit(log(hs),log(abs_error_b),1);
p(1) %%%% the slope is the estimate for the order of convergence


%%%% estimate of the order of convergence - centered f.d.
p=polyfit(log(hs),log(abs_error_c),1);
p(1) %%%% the slope is the estimate for the order of convergence



%%%% exercise 2

N = 10;
a = 0;
b = 1;

f=@(x)(4*cos(3*x));
exact=@(x)(1/9*(13 + 4*x*(-1 + cos(3)) - 4*cos(3*x)));

ua=exact(a);
ub=exact(b);

h=(b-a)/N;
xs=a+h*(1:(N-1));
fs=f(xs);

A = diag((-2/h^2)*ones(1,N-1)) + diag(1/(h^2)*ones(1,N-2),1) + diag(1/(h^2)*ones(1,N-2),-1);
bvec=fs';
bvec(1)=bvec(1)-ua/h^2;
bvec(end)=bvec(end)-ub/h^2;

us=A\bvec;


ts=linspace(a,b,100);


figure1 = figure(1);
axes1 = axes('Parent',figure1,'fontsize',16);
hold(axes1,'on');
plot(ts,exact(ts),'DisplayName','exact solution','LineWidth',2);
plot1 = plot(xs,us,'DisplayName','appr. solution','Marker','.','LineStyle','none','MarkerSize',17);
box(axes1,'on');
hold(axes1,'off');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.157589285714286 0.830952555934588 0.174107142857143 0.0630952380952381]);
xlim([a,b])



figure1 = figure(2);
axes1 = axes('Parent',figure1,'fontsize',16);
hold(axes1,'on');
plot1 = plot(xs,us-exact(xs'),'DisplayName','error','Marker','.','LineStyle','none','MarkerSize',17);
box(axes1,'on');
hold(axes1,'off');
legend1 = legend(axes1,'show');
set(legend1,...
    'Position',[0.157589285714286 0.830952555934588 0.174107142857143 0.0630952380952381]);
xlim([a,b])






%%%% part 2

NS=3:100;
er=[];
hs=[];
for N =NS
h=(b-a)/N;
hs=[hs,h];
xs=a+h*(1:(N-1));
fs=f(xs);
A = diag((-2/h^2)*ones(1,N-1)) + diag(1/(h^2)*ones(1,N-2),1) + diag(1/(h^2)*ones(1,N-2),-1);
bvec=fs';
bvec(1)=bvec(1)-ua/h^2;
bvec(end)=bvec(end)-ub/h^2;
us=A\bvec;
ex=exact(xs);
mx=max(abs(us(:)-ex(:)));
er=[er,mx];
end




p=polyfit(log(hs),log(er),1);
p(1) %%%% the slope is the estimate for the order of convergence


