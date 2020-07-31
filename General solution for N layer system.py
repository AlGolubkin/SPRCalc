from sympy import *
import cmath
import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

PI = cmath.pi
i = complex(0, 1)
tetha_SPR = math.radians(50)            # degrees
calc_wl = 400e-9                        # metres


n = [1.00027, 1.54, 0.05, 2.756, 1.786] # Лист показателей преломления:  0 - AIR
                                                                       # 1 - GLASS (~BK7)
                                                                       # 2 - SILVER (Ag) (n_Mi)
                                                                       # 3 - Compensator (Fe2O3)
                                                                       # 4 - Decompensator (Al2O3)
eps = [8.85e-12, 2.247, -15.243, 14.2, 0] # Лист диэлектрических проницаемостей:  0 - AIR
                                                                       # 1 - GLASS (~BK7)
                                                                       # 2 - SILVER (Ag) (n_Mi)
                                                                       # 3 - Compensator (Fe2O3)
                                                                       # 4 - Decompensator (Al2O3)
miu = [1, 1, 1, 1, 1] # Лист магнитных проницаемостей:  0 - AIR
                                                                        # 1 - GLASS (~BK7)
                                                                        # 2 - SILVER (Ag) (n_Mi)
                                                                        # 3 - Compensator (Fe2O3)
                                                                        # 4 - Decompensator (Al2O3)
d = [0, 5e-2, 50e-9, 50e-9, 51e-9]      # Толщины слоев:                  0 - AIR
                                                                        # 1 - GLASS (~BK7)
                                                                        # 2 - SILVER (Ag) (n_Mi)
                                                                        # 3 - Compensator (Fe2O3)
                                                                        # 4 - Decompensator (Al2O3)
a, b, step = 45e-9, 60e-9, 5e-10
d3_mas = np.linspace(a, b, (b-a)/step)

teta_res, n1, n2, n3, n4, n5, eps1, eps2, eps3, eps4, eps5, wl, d1, d2, d3, d4, d5 = symbols('teta_res n1 n2 n3 n4 n5 eps1 eps2 eps3 eps4 eps5 wl d1 d2 d3 d4 d5')

q1 = ((eps1-(n1**2)*(sin(teta_res)**2))**(1/2))/eps1
q2 = ((eps2-(n1**2)*(sin(teta_res)**2))**(1/2))/eps2
q3 = ((eps3-(n1**2)*(sin(teta_res)**2))**(1/2))/eps3
q4 = ((eps4-(n1**2)*(sin(teta_res)**2))**(1/2))/eps4
q5 = ((eps5-(n1**2)*(sin(teta_res)**2))**(1/2))/eps5

beta1 = (2*PI*d1/wl)*sqrt(eps1 - (n1**2)*sin(teta_res)**2)
beta2 = (2*PI*d2/wl)*sqrt(eps2 - (n1**2)*sin(teta_res)**2)
beta3 = (2*PI*d3/wl)*sqrt(eps3 - (n1**2)*sin(teta_res)**2)
beta4 = (2*PI*d4/wl)*sqrt(eps4 - (n1**2)*sin(teta_res)**2)
beta5 = (2*PI*d5/wl)*sqrt(eps5 - (n1**2)*sin(teta_res)**2)

M2 = Matrix([[cos(beta2), -i*sin(beta2/q2)], [-i*q2*sin(beta2), cos(beta2)]])
M3 = Matrix([[cos(beta3), -i*sin(beta3/q3)], [-i*q3*sin(beta3), cos(beta3)]])
M4 = Matrix([[cos(beta4), -i*sin(beta4/q4)], [-i*q4*sin(beta4), cos(beta4)]])

M = M2 * M3 * M4

M11 = M.row(0)[0]
M12 = M.row(0)[1]
M21 = M.row(1)[0]
M22 = M.row(1)[1] # Была ошибка: обращение к тому же элементу что и предыдущий : M21 а не M22

# Пусть знаменатель reflectivity coeffitient = 0

res = Eq((M11 + M12*q5)*q1 - (M21 + M22*q5), 0)
T = (M11 + M12*q5)*q1 - (M21 + M22*q5) # / (M11 + M12*q5)*q1 - (M21 + M22*q5)


# Примем что магнитная проницаемость во всех веществах равна 1, тогда eps = n**2
res_subs = []
for i in d3_mas:
    res_subs.append(abs(T.subs([(n1, n[1]), (eps1, n[1]**2), (eps2, eps[2]), (eps3, n[3]**2), (eps4, n[4]**2), (eps5, n[0]**2), (d2, d[2]), (d4, d[4]), (teta_res, tetha_SPR), (wl, calc_wl), (d3, i)])))
print(res_subs)

data = pd.DataFrame({'d3': d3_mas, 'res': res_subs})
r1 = data['res'].min()
print('min = ', data['d3'][data['res']==r1])

plt.plot(d3_mas, res_subs, label='Solutions')
plt.legend()
plt.show()

