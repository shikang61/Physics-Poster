import matplotlib.pyplot as plt
from scipy import constants
from scipy.optimize import curve_fit
from scipy import stats
import numpy as np


with open("data.txt", 'r') as f:
    name = []
    a_arc_sec = []
    a_arc_sec_error = []
    T_year = []
    T_year_error = []
    for L in f:
        L.strip()
        line = [L for L in L.split(" ") if L != ""]
        name.append(line[0])
        a_arc_sec.append(line[1])
        a_arc_sec_error.append(line[2])
        T_year.append(line[13])
        T_year_error.append(line[14])


bad_index = name.index('S111')
del name[bad_index]
del a_arc_sec[bad_index]
del a_arc_sec_error[bad_index]
del T_year[bad_index]
del T_year_error[bad_index]


def map_to_num(a):
    return list(map(float, a))

T_second = np.array(list(map(lambda x: x * constants.year, map_to_num(T_year))))
T_second_error = np.array(list(map(lambda x: x * constants.year, map_to_num(T_year_error))))
a_m = np.array(list(map(lambda x: x * constants.arcsec * 7940.0 * constants.parsec, map_to_num(a_arc_sec))))
a_m_error = np.array(list(map(lambda x: x * constants.arcsec * 7940.0 * constants.parsec, map_to_num(a_arc_sec_error)))) 


x = np.array([i for i in (T_second)])
y = np.array([i for i in (a_m)])
x_err = [2* (i/j) * k for i, j, k in zip(T_second_error, T_second, x)]
y_err = [3* (i/j) * k for i, j, k in zip(a_m_error, a_m, y)]

log_x = np.log(x)
log_y = np.log(y)
log_x_err = [i/j for i, j in zip(x_err, x)]
log_y_err = [i/j for i, j in zip(y_err, y)]

def func(x,a,b):
    return a*x + b


popt, pcov = curve_fit(func, log_x, log_y)
plt.scatter(log_x, log_y, marker='+', color='k', s=70)
plt.errorbar(log_x, log_y, yerr=log_y_err, xerr=log_x_err, fmt='none', 
    elinewidth=0.6, barsabove=True, capthick=0.4, capsize=4, color='k')
plt.plot(log_x, func(log_x, *popt), color='orange')
plt.xlabel("log(T / s)")
plt.ylabel("log(a / m)")
plt.title("Plot of ln(a) against ln(T)")


plt.xlim(19.5, 27)
plt.ylim(32.2, 38.5)
plt.annotate(rf'$ln\,T$ = {round(popt[0], 3)} $ln\,a$ + {round(popt[1], 5)}, $R^2$ = 0.99', xy=(0.1, 0.88), xycoords='axes fraction')
M = 1/np.exp(popt[1]*(-3)) * 4 * constants.pi**2 / (constants.gravitational_constant *1.98*10**36)
slope, intercept, r_value, p_value, std_err = stats.linregress(list(zip(log_x, log_y)))
print(std_err)
plt.annotate(fr'Mass of Black Hole = ${round((M), 3)} \pm {round(M*std_err, 4)}$ million $M_\odot$ ', xy=(0.1, 0.83), xycoords='axes fraction')

plt.show()
