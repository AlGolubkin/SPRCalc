from sympy import *
import cmath
import math
init_printing(pretty_print=True)

PI = cmath.pi
i = complex(0, 1)
tetha_SPR = math.radians(50)            # degrees
calc_wl = 400e-9                        # metres
n = [1.00027, 1.54, 0.05, 2.756, 1.786] # Лист показателей преломления:  0 - AIR
                                                                       # 1 - GLASS (~BK7)
                                                                       # 2 - SILVER (Ag) (n_Mi)
                                                                       # 3 - Compensator (Fe2O3)
                                                                       # 4 - Decompensator (Al2O3)
eps = [0, 2.247, -15.243, 14.2, 0] # Лист диэлектрических проницаемостей:  0 - AIR
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
# Найти произведение матриц Мi
x, y, beta_2, beta_3, beta_4, q2, q3, q4 = symbols('x y beta_2 beta_3 beta_4 q2 q3 q4')

M2 = Matrix([[cos(beta_2), -i*sin(beta_2/q2)],[-i*q2*sin(beta_2), cos(beta_2)]])
M3 = Matrix([[cos(beta_3), -i*sin(beta_3/q3)],[-i*q3*sin(beta_3), cos(beta_3)]])
M4 = Matrix([[cos(beta_4), -i*sin(beta_4/q4)],[-i*q2*sin(beta_4), cos(beta_4)]])

M = M2 * M3 * M4
M11 = M.row(0)[0]
M12 = M.row(0)[1]
M21 = M.row(1)[0]
M22 = M.row(1)[0]

M5 = M.det()

print('M11 = ' + str(M11))
print('M5 = ' + str(M5))
#print(simplify(M11))
#print(trigsimp(M11))
#print(solve(M5, beta_2))
#print(type(M))
#print(type(M11))
#print(M11.subs(beta_2, 2))

# Тест подстановки

# Попробуем найти толщину компенсирующего слоя из предыдущей статьи
# Попытка найти толщину из предыдущей статьи
#---------------------------
n_a, d_3, d_4, wl, n_Mi, n_3, n_4, n_1, eps_3, eps_4, mu_3, mu_4 = symbols('n_a d_3 d_4 wl n_Mi n_3 n_4 n_1 eps_3 eps_4 mu_3 mu_4')
n_eff = (n_Mi*(n_1*cmath.sin(tetha_SPR))/(n_Mi-(n_1*cmath.sin(tetha_SPR))**2))**(1/2)

n_b = Eq(((2*PI*d_3)/wl)*(((-n_Mi*(n_4**2))**(3/2))/(n_Mi-(n_4**2)))*(((n_3**2)-(n_4**2))/(n_3**2))+n_4,n_a)
n_b1 = n_b.subs([(n_a, n_eff), (n_3, 2.756), (n_4, 1.786), (n_Mi, -15.243), (wl, 400e-9), (n_1, 1.54)])
print(n_b1)
print('d3 = ',solve(n_b1, d_3))
print('sin(TETHA_spr) = ', math.sin(tetha_SPR))
#---------------------------
# Попробуем идею о равных бэта
eq1 = Eq((2*PI*d_3/wl)*((n_3**2 - (n_1**2)*(sin(tetha_SPR)**2))**(1/2)),(2*PI*d_4/wl)*((n_4**2 - (n_1**2)*(sin(tetha_SPR)**2))**(1/2)))
subs_eq1 = eq1.subs([(n_3, 2.756), (n_4, 1.786), (n_1, 1.54), (d_4, 50e-9)])
print('B3 = B4: ', subs_eq1)
eq1_1 = solve(subs_eq1, d_3)
print('d_3 = ' + str(eq1_1))
# НЕ РАБОТАЕТ ^
#---------------------------
# Попробуем сравнить beta_3 и beta_4
betha_3 = (2*PI*d_3/wl)*((n_3**2)-(n_1**2)*(sin(tetha_SPR)**2))
betha_4 = (2*PI*d_4/wl)*((n_4**2)-(n_1**2)*(sin(tetha_SPR)**2))

betha_3_sbs = betha_3.subs([(d_3, 50e-9), (wl, 400e-9), (n_3, 2.756), (n_1, 1.54)])
betha_4_sbs = betha_4.subs([(d_4, 51e-9), (wl, 400e-9), (n_4, 1.786), (n_1, 1.54)])

print(betha_3_sbs)
print(betha_4_sbs)
# НЕ СРАБОТАЛО ^

# Попробуем проверить самый первый вариант: получится ли вообще с помощью этого метода вычислить угол ППР

bth_k = lambda k : (2*PI*d[k]/calc_wl)*(n[k]**2 - (n[1]**2)*(math.sin(tetha_SPR)**2))**(1/2) # Расчет бэты
q_k = lambda k : (eps[k]**2 - (n[1]*(math.sin(tetha_SPR)**2))**(1/2)) / eps[k] # Что за sk?
print('beta(' + str(3) + ') = ' + str(bth_k(3)))


# Попробуем вычисления из статьи: High-perfomance bimetallic film surface plasmon resonance
#                                   sensor based on film thickness optimisation

