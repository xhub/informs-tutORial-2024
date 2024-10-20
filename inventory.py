import gamspy as gp
import gamspy.math as gpm
import numpy as np
import sys

# START PROBLEM DATA DEFINITION

m = gp.Container(system_directory="/Library/Frameworks/GAMS.framework/Resources")
# m = gp.Container()
i = gp.Set(m, name="i", records=['ordfix','inv'])
tau = gp.Parameter(m, name="tau", domain=i, records=np.array([10, 1000]), description="target")
alpha = gp.Parameter(m, name="alpha", domain=i, records=np.array([5, 0.5]))

rho = gp.Parameter(m, name='rho',description='units sold per year',records=144)
beta = gp.Parameter(m, name='beta',description='fixed cost of an order',records=400)
gamma = gp.Parameter(m, name='gamma',description='inv cost per unit per year',records=60)

# END PROBLEM DATA DEFINITION

# START OPTIMIZATION PROBLEM DEFINITION
f = gp.Variable(m, name="f", domain=i)
x = gp.Variable(m, name="x")

x.lo[:] = 10

defcost = gp.Equation(m, name="defcost", domain=i)
defcost[i] = f[i] == (beta*rho/x).where[i.sameAs('ordfix')] + (gamma*x/2).where[i.sameAs('inv')]
# END OPTIMIZATION PROBLEM DEFINITION

def reset(sym):
    sym.l[:] = 0
    sym.m[:] = 0

def ReSHOPAnnotation(m, s):
    return m.addGamsCode("EmbeddedCode ReSHOP:\n" + s + "\nendEmbeddedCode")

biobjective = gp.Model(m, name="biobjective", equations=[defcost], problem="emp")

# Default solve (alpha = [5,0.5])
ReSHOPAnnotation(m,"""
deffn f(i) defcost(i)
main: min SUM(i, alpha(i)*f(i)) x
""",
)
reset(x)
biobjective.solve(output=sys.stdout, solver="reshop")

# Default solve (alpha = [1, 1])
ReSHOPAnnotation(m,"""
deffn f(i) defcost(i)
main: min SUM(i, f(i)) x
""",
)
reset(x)
biobjective.solve(output=sys.stdout, solver="reshop")


# Second solve (alpha = [1,1])
ReSHOPAnnotation(m,"""
deffn f(i) defcost(i)
h0: MP('smax', f(i))
main: min h0.valFn x
""",
)
reset(x)
biobjective.solve(output=sys.stdout, solver="reshop")

# Third solve (alpha = [5,0.5])

defcostO = gp.Equation(m, name="defcostO", domain=i)
defcostO[i] = f[i] == (beta*rho/x).where[i.sameAs('ordfix')] + (gamma*x/2).where[i.sameAs('inv')] - tau[i]

biobjectiveO = gp.Model(m, name="biobjectiveO", equations=[defcostO], problem="emp")

ReSHOPAnnotation(m,"""
deffn f(i) defcostO(i)
H(i): MP('plus',f(i))
root: min SUM(i, alpha(i)*H(i).valfn) x
""",
)
reset(x)
biobjectiveO.solve(output=sys.stdout, solver="reshop")
