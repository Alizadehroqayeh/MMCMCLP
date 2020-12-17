scalar starttime; starttime = jnow;

option optcr=0;
set i/1*30/
j/1*300/
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
3  1   46    0;
*4  1   41    61;
table cp(l,k)
   1   2   3
1  40  70  100
2  60  110  0
3  45  100  0;
*4  40  60  110;
display xb,yb,dis,cp,d,g,cap,sd,h;

variable
obj,obj2;
binary variable
y(i,l,k,t),z(i);
positive variable
x(i,j,l,t);
equation
objfun,c2(i,l,t),c4(i,l,t),c5(i,l,t),c3(j,l,t),c6(i,t),c7;

objfun..       obj=e=sum((i,j,l,t),g(i,j)*d(j,l,t)*x(i,j,l,t))-sum((i,l,k,t),h(i,l,k,t)*y(i,l,k,t));
c2(i,l,t)..    sum(k,y(i,l,k,t))=l=z(i);

c3(j,l,t)..    sum(i,x(i,j,l,t))=l=1;
c4(i,l,t)..    sum(j,d(j,l,t)*x(i,j,l,t))=l=sum(k,cp(l,k)*y(i,l,k,t));
c5(i,l,t)..    sum(j,d(j,l,t)*x(i,j,l,t))=g=sum(k,cpl(l,k)*y(i,l,k,t));

c6(i,t)..      sum((l,j),d(j,l,t)*x(i,j,l,t))=l=cap(i)*z(i);
c7..           sum(i,z(i))=l=p;

model sc2/objfun,c2,c3,c4,c5,c6,c7/;

******************************************************
solve sc2 max obj use mip;
scalar elapsed; elapsed = (jnow - starttime)*24*3600;
display elapsed;
display obj.l,z.l,y.l,x.l;
