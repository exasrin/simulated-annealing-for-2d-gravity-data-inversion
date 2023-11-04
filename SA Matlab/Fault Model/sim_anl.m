function [xx,ff,dff]=sim_anl(f,x0,l,u,kmax,q,TolFun,model);
x=x0; fx=f;
xo=x; fo=fx;
if nargin < 7, TolFun=1e-8; end
if nargin < 6, q=1; end
if nargin < 5, kmax=500; end %max iteration number
xx=[];ff=[];dff=[];
for k=1:kmax
    Ti=1000*((0.99)^k);
    mu=10^(100*((k/kmax)^q));
    dx=mu_inv(2*rand(size(x))-1,mu).*(u-l);
    x1=x+dx;
    x1=(x1<l).*l+(l<= x1).*(x1<=u).*x1+(u< x1).*u;
    [gm]=forward_gravity(x1,model); %Fungsi Forward
    fx1=gm;
    Nf=length(fx1);
    df=1/Nf*(sum(sqrt((fx1-fx).^2))); 
    if df<0|rand<exp(Ti*df/(abs(fx)+eps)/TolFun);
        x=x1; fx=fx1;end
    if fx<fo, xo=x; fo=fx1;end
    fx11=fx1';
    f11=f';
    dfo=1/Nf*(sum(sqrt((fx11-f11).^2))); % Seleisih Obs dan Model Call
    dff=[dff dfo];
    xx=[xx x1]; % parameter solution of SA
    ff=[ff fx1];
end
