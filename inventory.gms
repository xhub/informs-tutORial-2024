SCALAR rho    'units sold per year'        /144/
       beta   'fixed cost of an order'     /400/
       gamma  'inv cost per unit per year' /60/;

* target
set i / ordfix, inv /;
parameter tau(i) target /ordfix 10, inv 1000 /;
parameter alpha(i) /ordfix 5, inv 0.5 /;

VARIABLES x, f(i);

EQUATIONS defcost(i);

defcost(i)..
  f(i) =e= (beta*rho/x)$sameas(i,'ordfix') + (gamma*x/2)$sameas(i,'inv');

model biobjective / defcost /;
x.lo = 10;

* This is the first variant in Subsection 3.1 alpha = [5, 0.5]
EmbeddedCode ReSHOP:
  deffn f(i) defcost(i)
  main: min SUM(i, alpha(i)*f(i)) x
endEmbeddedCode
 
option clear=x;
x.lo = 10;
solve biobjective using emp;

* This is the second variant in Subsection 3.1 alpha = [1 1]

* min-max solve
embeddedCode ReSHOP:
  deffn f(i) defcost(i)
  main: min h0.valFn x
  h0: MP('smax', f(i))
endEmbeddedCode
 
option clear=x;
x.lo = 10;
solve biobjective using emp;

* dual solve
embeddedCode ReSHOP:
  deffn f(i) defcost(i)
  main: min h0.dual().valFn x
  h0: MP('smax', f(i))
endEmbeddedCode
 
option clear=x;
x.lo = 10;
solve biobjective using emp;

* Equilibrium solve
embeddedCode ReSHOP:
  deffn f(i) defcost(i)
  main: min h0.objFn x
  h0: MP('smax', f(i))
  equil: Nash(main,h0)
endEmbeddedCode
 
option clear=x;
x.lo = 10;
solve biobjective using emp;

* third variant alpha = [5,0.5]
equation defcostO(i);

defcostO(i)..
  f(i) =e= (beta*rho/x)$sameas(i,'ordfix') + (gamma*x/2)$sameas(i,'inv') - tau(i);

* following solves as min-max
embeddedCode ReSHOP:
  deffn f(i) defcostO(i)
  H(i): MP('plus',f(i))
  root: min SUM(i, alpha(i)*H(i).valfn) x
endEmbeddedCode

model biobjectiveO / defcostO /;

option clear=x;
x.lo = 10;
solve biobjectiveO using emp;
