import math

import numpy as np
from pymoo.algorithms.moo.nsga2 import NSGA2
from pymoo.optimize import minimize
from pymoo.core.problem import ElementwiseProblem
# from pymoo.problems import get_problem
from pymoo.util.ref_dirs import get_reference_directions
from pymoo.visualization.scatter import Scatter


from hdh_calc import *

class NSGA_II_DCHDH(ElementwiseProblem):   # Class where the present problem is declared for NSGA-II, inheriting from one of pymoo's classes




     def __init__(self, mswl=1, mswu=100, twminl=298.15, twminu=313.15, twmaxl=328.15, twmaxu=343.15,
                 rhinl=75, rhinu=90, rhoul=75, rhouu=90, efdl=80, efdu=90, efhl=75, efhu=90):
        super().__init__(n_var=7,
                         n_obj=3,
                         n_ieq_constr=3,
                         xl=np.array([mswl, twminl, twmaxl, rhinl, rhoul, efdl, efhl]),
                         xu=np.array([mswu, twminu, twmaxu, rhinu, rhouu, efdu, efhu]))

        self.i = 1


     def _evaluate(self, x, out, *args, **kwargs):


        Caldef = intermediate_calculus()  # A desalinator type object is created to be computed

        # The object "Caldef" is now assigned the initial values

        Caldef.msw =(x[0])/100
        Caldef.tswid = x[1] - 273.15
        Caldef.tswih = x[2] - 273.15
        Caldef.RH_ih = x[3]/100
        Caldef.RH_oh = x[4]/100
        Caldef.Efd = x[5]/100
        Caldef.Efh = x[6]/100
        Caldef.Sal = 35
        Caldef.Pr_a = 101325



        n1=2.049         # Manufacturer's characteristics for the packed bed
        n2=-0.779
        n3=1
        n4=0.632

        u = 10          # Overall heat transfer coefficient (W/(m2*K))
        irad = 1000     # Solar radiation intensity (W/m2)
        coec = 0.45     # Coefficient of conversion of solar energy into heat (-)

        func = Caldef.hdhcalcul(0,0,15,15,5,u,n1,n2,n3,n4,irad,coec,0)



        f1 = func[0]
        f2 = func[1]
        f3 = func[2]
        g1 = 0.15 - ((func[4] + func[5]) / func[6])
        g2 = (func[8] / (100 * func[9] * func[10])) - 35
        #g3bat = 5 - ((func[11] * 3600) / func[7])
        g3 = func[3]


        print(self.i, round(f1, 2), round(f2, 2), round(f3, 2))  # round(vc, 2), round(fn, 2), round(ap, 2)
        self.i += 1

        out["F"] = [f1, f2, f3]
        out["G"] = [g1,g2,g3] #g3[g1, g2, g3]




# Algorithm execution
problem = NSGA_II_DCHDH()

algorithm = NSGA2(pop_size=20)

res = minimize(problem, algorithm, ("n_gen", 250), verbose=False, seed=1)   

print('Idx', ' f1=PWA      ', 'f2=1/RR      ', 'f3=EDA      ', 'MFRsw   ', 'tw,min ',  'tw,max ', 'RH_in  ', 'RH_ou  ', 'Ed  ', 'Eh  ')
for i in range(0, len(res.F)):
    msw = (res.X[i][0]) / 100
    twmin = (res.X[i][1]) - 273.15
    twmax = (res.X[i][2]) - 273.15
    rhin = (res.X[i][3]) / 100
    rhou = (res.X[i][4]) / 100
    efd = (res.X[i][5]) / 100
    efh = (res.X[i][6]) / 100

    print(i+1, ' ', res.F[i], round(res.X[i][0] / 100, 4), round(res.X[i][1] - 273.15, 2), round(res.X[i][2] - 273.15, 2)
          , round(res.X[i][3] / 100, 2), round(res.X[i][4] / 100, 2), round(res.X[i][5] / 100, 2), round(res.X[i][6] / 100, 2))

plot = Scatter()
plot.add(res.F, edgecolor="red", facecolor="none")
plot.show()

