set i/1*10/
j/1*100/
t/1*3/
l/1*3/
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
u2(i,t)
u3(i,l,t)
dx(i,j,l,t)
dd(i,j)
rhs(i,l,t)
srhs
yl(i,l,t)
gd(i),gdd(l,t),gdd2(l,t),gdd3(l,t),h3,h5,h7,ht,e1(i,j,t),dgy(i,j,l,t),gy2(i,j,l,t),gy(i,j,l,t);
xb(i)=uniform(0,1000);
xbb(j)=uniform(0,1000);
ybb(j)=uniform(0,1000);
yb(i)=uniform(0,1000);
dis(i,j)=(power(xb(i)-xbb(j),2)+power(yb(i)-ybb(j),2))**(1/2);
d(j,l,t)=uniform(1,5);
e(i,j)$(dis(i,j)<300)=1;
g(i,j)$(dis(i,j)<300)=1;
g(i,j)$(dis(i,j)<400and dis(i,j)>300) =0.75;
g(i,j)$(dis(i,j)<500 and dis(i,j)>400) =0.5;
g(i,j)$(dis(i,j)<600 and dis(i,j)>500) =0.25;
g(i,j)$(dis(i,j)>600)=0;
h(i,l,k,t)=uniform(40,60);
p=3;
cap(i)=uniform(80,100);

table cpl(l,k)
   1   2     3
1  1   21    36
2  1   26    0
3  1   26    0;
*4  1   21    31;
table cp(l,k)
   1   2   3
1  20  35  50
2  30  55  0
3  25  50  0;
*4  20  30  55;

variable
obj,obj2,bound,objlower,upperbound;
binary variable
y(i,l,k,t),z(i);
positive variable
x(i,j,l,t),x2(i,j,l,t);
equation
objfun,objfun1,objfun2,objfun12,objfun3,c6(i,t),c2(i,l,k,t),c1(i,j,l,t),c7,c8(j,l,t)
,c3(i,l,t),c4(i,l,t),c5(i,l,t),c41(i,l,t),c51(j,l,t),c61(i,t),c31(i,l,t),c81(i,j,l,t),c11(i,j,l,t);

objfun..         bound=e=sum((i,j,l,t),g(i,j)*d(j,l,t)*x(i,j,l,t))+sum((j,l,t),u(j,l,t)*(1-sum(i,x(i,j,l,t))))
+sum((i,j,l,t),u2(i,t)*(-d(j,l,t)*x(i,j,l,t)))-sum((i,l,k,t),h(i,l,k,t)*y(i,l,k,t))+sum((i,t),u2(i,t)*cap(i)*z(i));

c1(i,j,l,t)..   x(i,j,l,t)=l=sum(k,y(i,l,k,t));
c3(i,l,t)..    sum(k,y(i,l,k,t))=l=z(i);
c4(i,l,t)..     sum(j,d(j,l,t)*x(i,j,l,t))=l=sum(k,cp(l,k)*y(i,l,k,t));
c5(i,l,t)..     sum(j,d(j,l,t)*x(i,j,l,t))=g=sum(k,cpl(l,k)*y(i,l,k,t));
c7..            sum(i,z(i))=l=p;


model DMC/objfun,c1,c4,c5,c3,c7/;

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
parameter gamma(j,l,t),gamma2(i,t),gamma3(i,l,t),dy(i,j,l,t),maxdy(j,l,t),maxdx(j,l,t),b42(i,l,t),b52(i,l,t),b8(j,l,t),max2dx(i,l,t),xsum,x2lower,ylower,rhs(i,l,t),srhs,xlower;
scalar norm;
scalar lowerbound/80 /;

parameter uprevious(j,l,t),bound2,uprevious2(i,t),uprevious3(i,l,t),xl2(i,j,l,t),lowerbound2,lowerbound3,mainlower2,mainlower,maxjdx(i,l,t),xl(i,j,l,t),sumx(i,l,t),cap2(i,l,t),capl2(i,l,t);
scalar deltau;
parameter results(iter,*),zlower,y4;

*
* initialize u with relaxed duals
*
u(j,l,t) =0;
u2(i,t) =uniform(0.01,0.09);
*
* an upperbound on L
*

scalar starttime; starttime = jnow;
loop(iter$continue,
*
* solve the lagrangian dual problem
*

solve DMC max bound use mip;
if (bound.l < bestbound,
bestbound =bound.l;
results(iter,'upperbound') = bestbound;
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
gamma2(i,t) =sum((l,j),d(j,l,t)*x.l(i,j,l,t))-cap(i)*z.l(i) ;
norm = sum((j,l,t),sqr(gamma(j,l,t)))+ sum((i,l,t),sqr(gamma2(i,t)));
stepsize = theta*(bestbound-lowerbound)/norm;
*
* update duals u
*
uprevious(j,l,t) = u(j,l,t);
u(j,l,t) = max(0, u(j,l,t)+stepsize*gamma(j,l,t));
uprevious2(i,t) = u2(i,t);
u2(i,t) = max(0, u2(i,t)+stepsize*gamma2(i,t));
*
* converged ?
*

*);
if( (bestbound-lowerbound)/bestbound <0.01,
display "Converged";
continue = 0;
);
*if((results(iter,'upperbound')=results(iter-1,'upperbound')) and (results(iter-1,'upperbound')=results(iter-2,'upperbound')),
*display "converged";
*continue =0;
*);

);
scalar elapsed; elapsed = (jnow - starttime)*24*3600;
display elapsed;

display results,bound.l,x.l,z.l,y.l,u,u2;

