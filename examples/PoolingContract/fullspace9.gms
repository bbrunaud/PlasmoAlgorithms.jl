
SETS
i	"feeds"	/1*5/
j "products" /1*3/
l "pools" /1*4/
k "qualities" /1*2/
ww "base scenarios" /1*3/

;
Sets Tx(i, l), Ty(l, j), Tz(i, j);

Tx(i, l) = yes; Ty(l, j) = yes; Tz(i, j) = no;


PARAMETERS
AU(i)
/1 300, 2 250, 3 250, 4 200, 5 300 /,
AL(i)
/1 0, 2 0, 3 0, 4 0, 5 0 /,
SL(l)
/1 0, 2 0, 3 0, 4 0/,
SU(l)
/1 400, 2 400, 3 350, 4 500/,
base_d(j)
/1 5.7, 2  6.2, 3  6.8/,
base_DU(j)
/1 229, 2 173, 3 284/,
c_fixed_pool(l)
/1 310, 2 470, 3 380, 4 510/,
base_c_variable_pool(l)
/ 1 1.1, 2 0.9, 3 1.05, 4 0.8/,
c_variable_pool(l)
/1 1.1, 2 0.9, 3 1.05, 4 0.8/,
c_fixed_inlt(i)
/ 1 260, 2 70, 3 150, 4 190, 5 110/,
c_variable_inlt(i) 
/1 0.5, 2 0.8, 3 0.6, 4 0.55, 5 0.7/,
base_c_variable_inlt(i)
/1 0.5, 2 0.8, 3 0.6, 4 0.55, 5 0.7/
;

Table CC(i, k)		## Feed Concentration
      1     2
   1  0.13  0.87
   2  0.89  0.11 
   3  0.69  0.31
   4  0.28  0.72
   5  0.35  0.65;

Table PU(j,k)		## Maximum Allowable Product Concentration
     1     2
  1  0.56  0.44
  2  0.30  0.70
  3  0.41  0.59;

Table PL(j, k)	        ## Minimum Allowable Product Concentration
     1     2
  1  0.56  0.44
  2  0.30  0.70
  3  0.41  0.59;


Table c_x(i, l)		## Cost Parameter
	1	2	3	4
1	6.2     9.4     7.6     10.2
2	1.67    2.53    2.05    2.75
3	3.58    5.42    4.39    5.89
4	4.53    6.87    5.55    7.45
5	2.62    3.98    3.22    4.32;

PARAMETERS
base_prob(ww)
/1 0.3,2 0.4, 3 0.3/,
sigma_d(i),
sigma_b(i),
base_psi_f(i)
/1 0.5, 2 0.5, 3 0.5, 4 0.5, 5 0.5/,
base_psi_d1(i)
/1 0.55, 2 0.55, 3 0.55, 4 0.55, 5 0.55/,
base_psi_d2(i)
/1 0.4, 2 0.4, 3 0.4, 4 0.4, 5 0.4/,
base_psi_b1(i)
/1 0.55, 2 0.55, 3 0.55, 4 0.55, 5 0.55/,
base_psi_b2(i)
/1 0.48, 2 0.48, 3 0.48, 4 0.48, 5 0.48/,
ratio(ww)
/1 0.7, 2 1.0, 3 1.3/
;
sigma_b(i) = AU(i) / 2;
sigma_d(i) = AU(i) /3 * 2;
alias (ww, ww1, ww2, ww3);
sets 
w "scenarios" /1*9/;
Parameter
num, prob(w),
psi_f(i,w),
psi_d1(i,w),
psi_d2(i,w),
psi_b1(i,w),
psi_b2(i,w),
d(j, w),
DU(j,w);
alias(w, w1, w2);
loop(ww1,
  loop(ww2,
	num = (ord(ww1)-1)*3 + ord(ww2);
	prob(w2)$(ord(w2)=num) = base_prob(ww1) * base_prob(ww2); 
	DU(j, w2)$(ord(w2)=num) = base_DU(j) * ratio(ww1);
	d(j,w2)$(ord(w2)=num) = base_d(j);
	psi_f(i,w2)$(ord(w2)=num) = base_psi_f(i) * ratio(ww2);
	psi_d2(i,w2)$(ord(w2)=num) = base_psi_d2(i) * ratio(ww2);
	psi_d1(i,w2)$(ord(w2)=num) = base_psi_d1(i) * ratio(ww2) ;
	psi_b1(i,w2)$(ord(w2)=num) = base_psi_b1(i) * ratio(ww2) ; 
	psi_b2(i,w2)$(ord(w2)=num) = base_psi_b2(i) 	* ratio(ww2);
);
	);
POSITIVE VARIABLES
S(l)
A(i)
y(l, j, w)
z(i, j, w)
q(i, l, w)
CT(i, w)
CTf(i, w)
CTb(i, w)
CTd(i, w)
Bf(i, w)
Bd(i, w)
Bd1(i, w)
Bd2(i, w)
Bd11(i, w)
Bd12(i, w)
Bb(i, w)
Bb1(i, w)
Bb2(i, w);

Binary VARIABLES
uf(i, w)
ub(i, w)
ud(i, w)
ub1(i, w)
ub2(i, w)
ud1(i, w)
ud2(i, w)
gamma_intlt(i)
gamma_pool(l);

VARIABLES
obj ;
equations
f1(i),f2(i),f3(l),f4(l), e2(i,w),e3(l,w),e5(j,w),e6(l,w),e8(j,k,w),e9(j,k,w),
c1(i,w),c2(i,w), c3(i,w),c4(i,w),c5(i,w),c6(i,w),c7(i,w),c8(i,w),c9(i,w),c10(i,w),
c11(i,w),c12(i,w),c13(i,w),c14(i,w),c15(i,w),c16(i,w),c17(i,w),c18(i,w),c19(i,w),c20(i,w),
c21(i,w),c22(i,w),cobj;

f1(i) .. AL(i)*gamma_intlt(i)=l= A(i);
f2(i) .. A(i) =l= AU(i)*gamma_intlt(i);
f3(l) .. SL(l)* gamma_pool(l)=l=S(l);
f4(l) .. S(l)=l= SU(l)* gamma_pool(l);

	
e2(i, w) .. sum((l,j)$(Tx(i,l) and Ty(l,j)), q(i,l,w)*y(l,j,w))  + sum(j$Tz(i,j), z(i,j,w) )=l= A(i);
e3(l, w) .. sum(j$Ty(l,j), y(l,j,w) ) =l= S(l);

e5(j, w) .. sum(l$Ty(l,j), y(l,j,w)) + sum(i$Tz(i,j), z(i,j,w) )=l= DU(j,w);
e6(l, w) .. sum(i$Tx(i,l), q(i,l,w)) =e= 1 ;
	
e8(j, k, w) .. PL(j,k) * (sum(l$Ty(l,j), y(l,j,w)) + sum(i$Tz(i,j), z(i,j,w) )) =l= sum(i$Tz(i,j), CC(i,k) * z(i,j,w)) + sum((l,i)$(Ty(l,j) and Tx(i,l)), CC(i,k)*q(i,l,w)*y(l,j,w));
e9(j, k, w) .. PU(j,k) * (sum(l$Ty(l,j), y(l,j,w) ) + sum(i$Tz(i,j), z(i,j,w) )) =g= sum(i$Tz(i,j), CC(i,k) * z(i,j,w) ) + sum((l,i)$(Ty(l,j) and Tx(i,l)), CC(i,k)*q(i,l,w)*y(l,j,w) );

c1(i, w) .. sum((l,j)$(Tx(i,l) and Ty(l,j)), q(i,l,w)*y(l,j,w) )  + sum(j$Tz(i,j), z(i,j,w))=e= Bf(i,w)+Bd(i,w)+Bb(i,w);
c2(i, w) .. AL(i)*uf(i,w)=l=Bf(i,w);
c3(i, w) .. AU(i)*uf(i,w)=g=Bf(i,w);
c4(i, w) .. AL(i)*ud(i,w)=l=Bd(i,w);
c5(i, w) .. AU(i)*ud(i,w)=g=Bd(i,w);	
c6(i, w) .. AL(i)*ub(i,w)=l=Bb(i,w);
c7(i, w) .. AU(i)*ub(i,w)=g=Bb(i,w);
c8(i, w) .. uf(i,w)+ub(i,w)+ud(i,w)=l= gamma_intlt(i);
c9(i, w) .. CT(i,w)=e=CTf(i,w)+CTb(i,w)+CTd(i,w);

c10(i, w) .. CTf(i,w)=e=psi_f(i,w)*Bf(i,w);

c11(i, w) .. CTd(i,w)=e=psi_d1(i,w)*Bd1(i,w)+psi_d2(i,w)*Bd2(i,w);
c12(i, w) .. Bd(i,w)=e= Bd1(i,w)+Bd2(i,w);
c13(i, w) .. Bd1(i,w)=e= Bd11(i,w)+Bd12(i,w);
c14(i, w) .. Bd11(i,w)=l= sigma_d(i)*ud1(i,w);
c15(i, w) .. Bd12(i,w)=e=sigma_d(i)*ud2(i,w);
c16(i, w) .. Bd2(i,w)=l= AU(i)*ud2(i,w);


c17(i, w) .. CTb(i,w)=e= psi_b1(i,w)*Bb1(i,w)+psi_b2(i,w)*Bb2(i,w);
c18(i, w) .. Bb(i,w)=e= Bb1(i,w)+Bb2(i,w);
c19(i, w) .. Bb1(i,w)=l= sigma_b(i)*ub1(i,w);
c20(i, w) .. Bb2(i,w)=g= sigma_b(i)*ub2(i,w);
c21(i, w) .. Bb2(i,w)=l= AU(i)*ub2(i,w);
c22(i, w) .. ub1(i,w)+ub2(i,w)=e=ub(i,w);

cobj ..  obj=e= sum(i, c_fixed_inlt(i) * gamma_intlt(i)+c_variable_inlt(i)*A(i) ) + sum(l, c_fixed_pool(l) * gamma_pool(l) + c_variable_pool(l)*S(l)) + sum(w, prob(w)*(sum(i, CT(i,w) ) -sum(j, d(j,w)*(sum(l$Ty(l,j), y(l,j,w) ) + sum(i$Tz(i,j), z(i,j,w) )) ) ));


model pooling /all/;
OPTION LIMROW = 0;
OPTION LIMCOL = 0;
OPTION OPTCA  = 1E-09;
OPTION OPTCR  = 1E-03;
OPTION RESLIM = 1E+04;
OPTION ITERLIM = 1E+09;

OPTION LP=CPLEX;
OPTION NLP=SNOPT;
OPTION MIP=CPLEX;
*OPTION MINLP=bonmin;
*OPTION MINLP=ANTIGONE;
*OPTION MINLP=BARON;
*OPTION MINLP=COUENNE;
OPTION MINLP=SCIP;
solve pooling using MINLP minmizing obj;












