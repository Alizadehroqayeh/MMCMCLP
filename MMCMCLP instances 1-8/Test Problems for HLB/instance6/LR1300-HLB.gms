set i/1*30/
j/1*300/
t/1*2/
l/1*4/
k/1*3/;
parameter
c(i)
d(j,l,t)
xb(i)
xbb(j)
yb(i)
ybb(j)
dis(i,j)
e(i,j)
g(i,j)
cap(i)
ee(i,j,t)
eee(i,j,t)
h(i,l,k,t)
q(l)
oj
p
e1(i,j,t)
u(j,l,t)
u2(i,l,t)
u3(i,l,t)
dx(i,j,l,t)
dd(i,j)
rhs(i,l,t)
srhs
yl(i,l,t)
sd(t);
xb(i)=uniform(0,1000);
xbb(j)=uniform(0,1000);
ybb(j)=uniform(0,1000);
yb(i)=uniform(0,1000);
dis(i,j)=(power(xb(i)-xbb(j),2)+power(yb(i)-ybb(j),2))**(1/2);
d(j,l,t)=uniform(1,5);
sd(t)=sum((j,l),d(j,l,t));
e(i,j)$(dis(i,j)<300)=1;
g(i,j)$(dis(i,j)<300)=1;
g(i,j)$(dis(i,j)<400and dis(i,j)>300) =0.75;
g(i,j)$(dis(i,j)<500 and dis(i,j)>400) =0.5;
g(i,j)$(dis(i,j)<600 and dis(i,j)>500) =0.25;
g(i,j)$(dis(i,j)>600)=0;
h(i,l,k,t)=uniform(40,60);
p=12;
cap(i)=uniform(150,180);
table cpl(l,k)
   1   2     3
1  1   41    71
2  1   61    0
3  1   46    0
4  1   41    61;
table cp(l,k)
   1   2   3
1  40  70  100
2  60  110  0
3  45  100  0
4  40  60  110;
variable
obj,obj2,bound,objlower,upperbound;
binary variable
y(i,l,k,t),z(i);
positive variable
x(i,j,l,t),x2(i,j,l,t);
equation
objfun,c6(i,t),c2(i,l,t),c1(i,j,l,t),c7;

objfun..         bound=e=-sum((i,l,k,t),h(i,l,k,t)*y(i,l,k,t))+sum((i,j,l,t),g(i,j)*d(j,l,t)*x(i,j,l,t))+sum((j,l,t),u(j,l,t)*(1-sum(i,x(i,j,l,t))))
+sum((i,l,t),u2(i,l,t)*(-sum(j,d(j,l,t)*x(i,j,l,t))+sum(k,cp(l,k)*y(i,l,k,t))))+sum((i,l,t),u3(i,l,t)*(sum(j,d(j,l,t)*x(i,j,l,t))-sum(k,cpl(l,k)*y(i,l,k,t))));


c1(i,j,l,t)..  x(i,j,l,t)=l=sum(k,y(i,l,k,t));
c2(i,l,t)..    sum(k,y(i,l,k,t))=l=z(i);
c6(i,t)..      sum((j,l),d(j,l,t)*x(i,j,l,t))=l=cap(i)*z(i);
c7..           sum(i,z(i))=l=p;
model ldual/objfun,c2,c7,c6,c1/;

***************************************
*---------------------------------------------------------------------
* subgradient iterations
*---------------------------------------------------------------------
set iter /iter1*iter80/;
scalar continue /1/;
scalar continue2 /1/;
parameter stepsize;
scalar theta /2/;
scalar noimprovement /0/;
scalar bestbound /+INF/;
parameter gamma(j,l,t),gamma2(i,l,t),gamma3(i,l,t),dy(i,j,l,t),maxdy(j,l,t),maxdx(j,l,t),b42(i,l,t),b52(i,l,t),b8(j,l,t),max2dx(i,l,t),xsum,x2lower,ylower,rhs(i,l,t),srhs,xlower;
scalar norm;

scalar lowerbound/1819/;
parameter uprevious(j,l,t),bound2,uprevious2(i,l,t),uprevious3(i,l,t),xl2(i,j,l,t),lowerbound2,lowerbound3,mainlower2,mainlower,maxjdx(i,l,t),xl(i,j,l,t),sumx(i,l,t),cap2(i,l,t),capl2(i,l,t);
scalar deltau;
parameter results(iter,*),zlower,y4;

*
* initialize u with relaxed duals
*
u(j,l,t) =0;
u2(i,l,t) = 0;
u3(i,l,t) = 0;
*
* an upperbound on L
*
scalar starttime; starttime = jnow;
loop(iter$continue,
*
* solve the lagrangian dual problem
*
*option optcr=0;
solve ldual max bound using mip;
if (bound.l < bestbound,
bestbound = bound.l;
results(iter,'upperbound') = bestbound;

*solve DMC max obj use lp;
*if (obj.l > lowerbound,
*lowerbound =obj.l;
*results(iter,'lowerbound') = lowerbound;

*);
noimprovement = 0;
else
noimprovement = noimprovement + 1;
if (noimprovement > 1,
theta = theta*0.4;
noimprovement = 0;
);
);
*
* calculate step size
*
gamma(j,l,t) = sum(i,x.l(i,j,l,t))-1;
gamma2(i,l,t) =sum(j,d(j,l,t)*x.l(i,j,l,t))-sum(k,cp(l,k)*y.l(i,l,k,t)) ;
gamma3(i,l,t) =-sum(j,d(j,l,t)*x.l(i,j,l,t))+sum(k,cpl(l,k)*y.l(i,l,k,t));
norm = sum((j,l,t),sqr(gamma(j,l,t)))+ sum((i,l,t),sqr(gamma2(i,l,t)))+ sum((i,l,t),sqr(gamma3(i,l,t)));
stepsize = theta*(bestbound-lowerbound)/norm;
*
* update duals u
*
uprevious(j,l,t) = u(j,l,t);
u(j,l,t) = max(0, u(j,l,t)+stepsize*gamma(j,l,t));
uprevious2(i,l,t) = u2(i,l,t);
u2(i,l,t) = max(0, u2(i,l,t)+stepsize*gamma2(i,l,t));
uprevious3(i,l,t) = u3(i,l,t);
u3(i,l,t) = max(0, u3(i,l,t)+stepsize*gamma3(i,l,t));
*
* converged ?
*
deltau = smax((j,l,t),abs(uprevious(j,l,t)-u(j,l,t)));

if( (bestbound-lowerbound)/bestbound <0.1,
display "Converged";
continue = 0;
);
);
scalar elapsed; elapsed = (jnow - starttime)*24*3600;
display elapsed;

display results,x.l,z.l,y.l;

