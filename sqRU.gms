$ontext

Simple example to illustrate how to model the CVaR as a support function

author: Olivier Huber <oli.huber@gmail.com>

$offtext

* START PROBLEM DATA DEFINITION

$if not set nSamples $set nSamples 950

SETS i 'realizations set' /1*%nSamples%/
     j 'dimension of x'   /1*2/
     k 'number of distributions' /1*4/;




PARAMETERS xi(k,i) 'normal distributions',
           mean(k)   /1 250, 2  125, 3 2500, 4 40000/,
           stddev(k) /1  75, 2 62.5, 3  500, 4  4000/;


embeddedCode Python:

import numpy as np
np.random.seed(0)
k = list(gams.get("k"))
i = list(gams.get("i"))

si = len(i)
sk = len(k)

means   = dict(gams.get("mean"))
stddevs = dict(gams.get("stddev"))


xi = []
for kk in k:
   dvec = [kk]*si
   xi.extend(zip(dvec, i, np.random.normal(loc=means[kk], scale=stddevs[kk], size=si)))

# import IPython
# IPython.embed(colors="neutral")
gams.set("xi", xi)

endEmbeddedCode xi


SCALAR tail;
tail = 1-.05;

parameter p_a(i);
p_a(i) = 1/(card(i)*(1-tail));

* END PROBLEM DATA DEFINITION

* START OPTIMIZATION PROBLEM DEFINITION

VARIABLES phi(i), x(j), gamma, obj;

x.lo('1') = 0.1;
x.up('1') = 0.2;
x.lo('2') = 0.1;
x.up('2') = 0.6;

x.l(j) = x.lo(j);

EQUATIONS defphi(i), defObj;

defphi(i).. phi(i) =E= 4 * xi('1', i) / (x('1')*x('2')*xi('4', i)) + 4*xi('2', i) * x('2') / (sqr(x('1')) * xi('4', i)) + sqr(xi('3', i)) / (sqr(x('1')) * sqr(xi('4', i))) - gamma;
 
defObj.. gamma =E= obj;

* END OPTIMIZATION PROBLEM DEFINITION

$macro reset(x) x.l(j) = 0.; x.m(j) = 0.;

model superquantile_obj /defphi, defObj/;

EmbeddedCode ReSHOP:
deffn phi(i) defphi(i)
H(i): MP("plus", phi(i))
root: min obj + SUM(i, p_a(i)*H(i).valfn) x(j) gamma defObj
endEmbeddedCode

reset(x)
solve superquantile_obj using emp;


