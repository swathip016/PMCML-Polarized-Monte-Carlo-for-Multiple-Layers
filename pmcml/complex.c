/*2:*/
#line 9 "complex.w"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#include "complex.h"/*6:*//*5:*/
#line 116 "complex.w"
static void complex_error(char*s)/*:5*/
#line 120 "complex.w"
{
printf("%s\n",s);
exit(1);
}/*:6*//*8:*//*7:*/
#line 129 "complex.w"
struct complex cset(double a,double b)/*:7*/
#line 133 "complex.w"
{
struct complex c;
c.re=a;
c.im=b;
return c;
}/*:8*//*10:*//*9:*/
#line 144 "complex.w"
struct complex cpolarset(double r,double theta)/*:9*/
#line 148 "complex.w"
{
return cset(r*cos(theta),r*sin(theta));
}/*:10*//*12:*//*11:*/
#line 157 "complex.w"
double cabbs(struct complex z)/*:11*/
#line 161 "complex.w"
{
double x,y,temp;

x=fabs(z.re);
y=fabs(z.im);
/*printf("x=%f,y=%f",x,y);*/
if(x==0.0)return y;
if(y==0.0)return x;

if(x>y){
temp=y/x;
return(x*sqrt(1.0+temp*temp));

}

temp=x/y;
return(y*sqrt(1.0+temp*temp));
}/*:12*//*16:*//*15:*/
#line 190 "complex.w"
double carg(struct complex z)/*:15*/
#line 194 "complex.w"
{
return atan2(z.im,z.re);
}/*:16*//*18:*//*17:*/
#line 201 "complex.w"
double cnorm(struct complex z)/*:17*/
#line 205 "complex.w"
{
return(z.re*z.re+z.im*z.im);
}/*:18*//*20:*//*19:*/
#line 210 "complex.w"
struct complex csqrt(struct complex z)/*:19*/
#line 214 "complex.w"
{
double a,b;

if((z.re==0.0)&&(z.im==0.0))
return cset(0,0);

a=sqrt((fabs(z.re)+cabbs(z))*0.5);
if(z.re>=0)
b=z.im/(a+a);
else{
b=z.im<0? -a:a;
a=z.im/(b+b);
}

return cset(a,b);
}/*:20*//*22:*//*21:*/
#line 235 "complex.w"
struct complex csqr(struct complex z)/*:21*/
#line 239 "complex.w"
{
return cmul(z,z);
}/*:22*//*24:*//*23:*/
#line 246 "complex.w"
struct complex cinv(struct complex w)/*:23*/
#line 250 "complex.w"
{
double r,d;

if((w.re==0)&&(w.im==0))complex_error("Attempt to invert 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r=w.im/w.re;
d=1/(w.re+r*w.im);
return cset(d,-r*d);
}

r=w.re/w.im;
d=1/(w.im+r*w.re);
return cset(r*d,-d);
}/*:24*//*14:*//*13:*/
#line 181 "complex.w"
struct complex conj(struct complex z)/*:13*/
#line 185 "complex.w"
{
return cset(z.re,-z.im);
}/*:14*//*27:*//*26:*/
#line 271 "complex.w"
struct complex cadd(struct complex z,struct complex w)/*:26*/
#line 275 "complex.w"
{
struct complex c;

c.im=z.im+w.im;
c.re=z.re+w.re;
return c;
}/*:27*//*29:*//*28:*/
#line 286 "complex.w"
struct complex csub(struct complex z,struct complex w)/*:28*/
#line 290 "complex.w"
{
struct complex c;
c.im=z.im-w.im;
c.re=z.re-w.re;
return c;
}/*:29*//*31:*//*30:*/
#line 300 "complex.w"
struct complex cmul(struct complex z,struct complex w)/*:30*/
#line 304 "complex.w"
{
struct complex c;
c.re=z.re*w.re-z.im*w.im;
c.im=z.im*w.re+z.re*w.im;
return c;
}/*:31*//*33:*//*32:*/
#line 315 "complex.w"
struct complex cdiv(struct complex z,struct complex w)/*:32*/
#line 319 "complex.w"
{
struct complex c;
double r,denom;

if((w.re==0)&&(w.im==0))complex_error("Attempt to divide by 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r=w.im/w.re;
denom=w.re+r*w.im;
c.re=(z.re+r*z.im)/denom;
c.im=(z.im-r*z.re)/denom;
}else{
r=w.re/w.im;
denom=w.im+r*w.re;
c.re=(z.re*r+z.im)/denom;
c.im=(z.im*r-z.re)/denom;
}
return c;
}/*:33*//*35:*//*34:*/
#line 343 "complex.w"
double crdiv(struct complex z,struct complex w)/*:34*/
#line 347 "complex.w"
{
double r,c,denom;

if((w.re==0)&&(w.im==0))complex_error("Attempt to find real part with divisor 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r=w.im/w.re;
denom=w.re+r*w.im;
c=(z.re+r*z.im)/denom;
}else{
r=w.re/w.im;
denom=w.im+r*w.re;
c=(z.re*r+z.im)/denom;
}
return c;
}/*:35*//*37:*//*36:*/
#line 367 "complex.w"
double crmul(struct complex z,struct complex w)/*:36*/
#line 371 "complex.w"
{
return z.re*w.re-z.im*w.im;
}/*:37*//*42:*//*41:*/
#line 394 "complex.w"
struct complex csadd(double x,struct complex z)/*:41*/
#line 398 "complex.w"
{
struct complex c;
c.re=x+z.re;
c.im=z.im;
return c;
}/*:42*//*44:*//*43:*/
#line 409 "complex.w"
struct complex csdiv(double x,struct complex w)/*:43*/
#line 413 "complex.w"
{
struct complex c;
double r,factor;

if((w.re==0)&&(w.im==0))complex_error("Attempt to divide scalar by 0+0i");

if(fabs(w.re)>=fabs(w.im)){
r=w.im/w.re;
factor=x/(w.re+r*w.im);
c.re=factor;
c.im= -r*factor;
}else{
r=w.re/w.im;
factor=x/(w.im+r*w.re);
c.im= -factor;
c.re=r*factor;
}
return c;
}/*:44*//*40:*//*39:*/
#line 380 "complex.w"
struct complex csmul(double x,struct complex z)/*:39*/
#line 384 "complex.w"
{
struct complex c;
c.re=z.re*x;
c.im=z.im*x;
return c;
}/*:40*//*47:*//*46:*/
#line 438 "complex.w"
struct complex csin(struct complex z)/*:46*/
#line 442 "complex.w"
{
return cset(sin(z.re)*cosh(z.im),cos(z.re)*sinh(z.im));
}/*:47*//*49:*//*48:*/
#line 449 "complex.w"
struct complex ccos(struct complex z)/*:48*/
#line 453 "complex.w"
{
return cset(cos(z.re)*cosh(z.im),-(sin(z.re)*sinh(z.im)));
}/*:49*//*51:*//*50:*/
#line 480 "complex.w"
struct complex ctan(struct complex z)/*:50*/
#line 484 "complex.w"
{
double t,x,y;

if(z.im==0)return cset(tan(z.re),0.0);
if(z.im>DBL_MAX_10_EXP)return cset(0.0,1.0);
if(z.im< -DBL_MAX_10_EXP)return cset(0.0,-1.0);

x=2*z.re;
y=2*z.im;
t=cos(x)+cosh(y);
if(t==0)complex_error("Complex tangent is infinite");

return cset(sin(x)/t,sinh(y)/t);
}/*:51*//*53:*//*52:*/
#line 502 "complex.w"
struct complex casin(struct complex z)/*:52*/
#line 506 "complex.w"
{
struct complex x;
x=clog(cadd(cset(-z.im,z.re),csqrt(csub(cset(1,0),cmul(z,z)))));
return cset(x.im,-x.re);
}/*:53*//*55:*//*54:*/
#line 515 "complex.w"
struct complex cacos(struct complex z)/*:54*/
#line 519 "complex.w"
{
struct complex x;
x=clog(cadd(z,cmul(cset(0,1),csqrt(csub(cset(1,0),csqr(z))))));
return cset(x.im,-x.re);
}/*:55*//*57:*//*56:*/
#line 528 "complex.w"
struct complex catan(struct complex z)/*:56*/
#line 532 "complex.w"
{
struct complex x;
x=clog(cdiv(cset(z.re,1+z.im),cset(-z.re,1-z.im)));
return cset(-x.im/2,x.re/2);
}/*:57*//*62:*//*61:*/
#line 550 "complex.w"
struct complex csinh(struct complex z)/*:61*/
#line 554 "complex.w"
{
return cset(sinh(z.re)*cos(z.im),cosh(z.re)*sin(z.im));
}/*:62*//*60:*//*59:*/
#line 541 "complex.w"
struct complex ccosh(struct complex z)/*:59*/
#line 545 "complex.w"
{
return cset(cosh(z.re)*cos(z.im),sinh(z.re)*sin(z.im));
}/*:60*//*64:*//*63:*/
#line 559 "complex.w"
struct complex ctanh(struct complex z)/*:63*/
#line 563 "complex.w"
{
double x=2*z.re;
double y=2*z.im;
double t=1.0/(cosh(x)+cos(y));

return cset(t*sinh(x),t*sin(y));
}/*:64*//*66:*//*65:*/
#line 572 "complex.w"
struct complex catanh(struct complex z)/*:65*/
#line 576 "complex.w"
{
return catan(cset(-z.im,z.re));
}/*:66*//*68:*//*67:*/
#line 581 "complex.w"
struct complex casinh(struct complex z)/*:67*/
#line 585 "complex.w"
{
return casin(cset(-z.im,z.re));
}/*:68*//*71:*//*70:*/
#line 592 "complex.w"
struct complex cexp(struct complex z)/*:70*/
#line 596 "complex.w"
{
double x=exp(z.re);
return cset(x*cos(z.im),x*sin(z.im));
}/*:71*//*73:*//*72:*/
#line 602 "complex.w"
struct complex clog(struct complex z)/*:72*/
#line 606 "complex.w"
{
return cset(log(cabbs(z)),carg(z));
}/*:73*//*75:*//*74:*/
#line 611 "complex.w"
struct complex clog10(struct complex z)/*:74*/
#line 615 "complex.w"
{
return cset(0.2171472409516259*log(cnorm(z)),carg(z));
}/*:75*//*78:*//*77:*/
#line 624 "complex.w"
struct complex*new_carray(long size)/*:77*/
#line 628 "complex.w"
{
struct complex*a;

if(size<=0)complex_error("Non-positive complex array size chosen");

a=(struct complex*)malloc(size*sizeof(struct complex));

if(a==NULL)complex_error("Can't allocate complex array");
return a;
}/*:78*//*80:*//*79:*/
#line 640 "complex.w"
void free_carray(struct complex*a)/*:79*/
#line 644 "complex.w"
{
if(a!=NULL)free(a);
}/*:80*//*82:*//*81:*/
#line 652 "complex.w"
struct complex*copy_carray(struct complex*a,long size)/*:81*/
#line 656 "complex.w"
{
struct complex*b=NULL;

if(a==NULL)complex_error("Can't duplicate a NULL complex array");

b=new_carray(size);
if(b!=NULL)memcpy(b,a,size*sizeof(struct complex));
return b;
}/*:82*//*84:*//*83:*/
#line 669 "complex.w"
void set_carray(struct complex*a,long size,struct complex z)/*:83*/
#line 673 "complex.w"
{
long j;

if(a==NULL)complex_error("Can't operate on a NULL complex array");

for(j=0;j<size;j++)a[j]=z;
}/*:84*//*:2*/
