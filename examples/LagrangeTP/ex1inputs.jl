using JuMP
using Gurobi

oproducts = 1:3
#m
markets = 1:3
#s
sites = 1:3
#t
otime = 1:3
#nt

nt = 3
#horas disponibles de produccion (longitud del periodo de tiempo)
HT= 720

#table sut(s,i)
sut=[
  100    60      80
  70     100     50
  80     40      100 ]

#  Table prate(s,i) production rate kg per hr
prate =[
     .3       .9        .3
     .1       .5        .5
     .1       .1        .2]

#table alpha(s,i)  Production cost of product i in production site s during time period t -!= f(t)
α = [
             5.5	2.5	2.5
             3.5	4.5	5.5
             4.5	5.5	4.5]

#table fcast(m,i,t)
fcast = zeros(3,3,3)
fcast[:,1,:] = [
             70	120	17
             120	90	8
             50	80	8]
fcast[:,2,:] = [
             90	140	20
             100	80	7
             500	70	6]
fcast[:,3,:] = [
             170	120	13
             600	90	5
             100	80	15]


#table beta(m,i) sale price
#bb(m,i) sale price of product i in market m during time period t
β=[
             20	23	26
             25	28	32
             30	33	37]

#table gamma(m,s,i,t) shipment cost- != f(t)
γ = zeros(3,3,3)
γ[1,:,:]=[
3	3	3
4	4	2
1	1	3]
γ[2,:,:]=[
2	2	2
4	1	4
5	3	2]
γ[3,:,:]=[
5	5	5
2	3	2
1	1	5]

#table delta(s,i,t) inventory cost- != f(t)
δ=[0.1	0.15	0.1
  0.2	0.2	0.25
  0.05	0.15	0.1]

#table setupcost(s,i,t)- != f(t)
setupcost=[
             100	200	100
             400	300	200
             700	500	400]
