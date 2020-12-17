scalar starttime; starttime = jnow;
option optcr=0;
set i/1*10/
j/1*100/
t/1*2/
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
permit(i,l,t)
gd(i),gdd(l,t),gdd2(l,t),gdd3(l,t),h3,h5,h7,ht,e1(i,j,t);
permit(i,l,t)=no;
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
2  1   31    0
3  1   26    0;
*4  1   21    31;
table cp(l,k)
   1   2   3
1  20  35  50
2  30  55  0
3  25  50  0 ;
*4  20  30  55;
display cap;
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






