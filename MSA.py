#! /usr/bin/env python
# -*- coding: utf-8 -*-

########################################################################################
#Определение наклона и позиционного угла спиральных галактик методом Полторака-Фридмана#
########################################################################################

from math import pi, radians, degrees, isnan
import sys
from os import system

try:
    from pylab import *
except ImportError:
    print "Error (03): can not import name pylab"
    sys.exit(1)

try:
    from numpy import empty, array, hypot, arange, mean, std, cos, sin
except ImportError:
    print "Error (04): can not import name numpy"
    sys.exit(1)

def ismono(r):                          # функция, возвращающая 1, если последовательность монотонно возрастает (строго)
    for i in xrange(len(r)-1):          # и 0 в противном случае
        if r[i] >= r[i+1]:
            return 0
    return 1


def rotate(i, PA, x, y):                          # Функция, аффинным преобразованием моделирующая 
    x = array(x)
    y = array(y)
    sinPA, cosPA = sin(PA), cos(PA)               # различные ориентации галактики в пространстве
    x1 = x*cosPA + y*sinPA      # Поворот декартовой системы координат вокруг
    y1 = y*cosPA - x*sinPA      # начала на угол PA
    x1 /= cos(i)                # Растяжение компоненты x               
    return x1, y1   

nArms = int(sys.argv[2]) # number of arms

#system('inmidas -j "@@ MSA.prg %s %i"' % (sys.argv[1], nArms))

pfile = open('tab_center.dat').readlines()
cenx_p = float(pfile[-2].split()[1])
ceny_p = float(pfile[-2].split()[2])


pointsOnArms = [] # list with points of all arms
for arm in xrange(1, nArms + 1):
    fileArm = open('arm%i.dat'%(arm), 'rb')
    pointsOnArms.append([])
    for s in fileArm.readlines()[3:-1]:
        x = float(s.split()[1]) - cenx_p
        y = float(s.split()[2]) - ceny_p
        pointsOnArms[arm-1].append((x, y))
    if (len(pointsOnArms[arm-1]) < 3):
        print "Too few points on arm %i; exiting..." % (arm)
        sys.exit(1)
    fileArm.close()

true_i, true_PA = [], []
for arm in xrange(nArms):
    true_i.append([])
    true_PA.append([])
    currentArmX = [i[0] for i in pointsOnArms[arm]]
    currentArmY = [i[1] for i in pointsOnArms[arm]]
    for incl in arange(radians(1), radians(90), radians(0.5)): # пробегаем все возможные углы
        for PA in arange(0, radians(180), radians(0.5)): # с шагом пол градуса
            r = hypot(*rotate(incl, PA, currentArmX, currentArmY)) # расстояние точек рукава от центра галактики (после поворота)
            if ismono(r) == 1:
                true_i[arm].append(degrees(incl))
                true_PA[arm].append(degrees(PA))


isConnected = [True] * nArms
separators = zeros(nArms) + 90 # По умолчанию разделители находятся на 90 град.
for arm in xrange(nArms):
    true_PA_sort = sorted(true_PA[arm])
    for i in xrange(len(true_PA_sort)-1):
        if true_PA_sort[i+1] - true_PA_sort[i] > 1.9:
            isConnected[arm] = False
            separators[arm] = (true_PA_sort[i] + true_PA_sort[i+1]) / 2.0
            break
            
#print isConnected
#print separators

true_PA_sep = true_PA[:]
# Если хотябы одно из множеств разомкнуто
if not all(isConnected):
    # Сдвигаем все множества в диапазон -90:90
    for arm in xrange(nArms):
        for i in xrange(len(true_PA[arm])):
            if true_PA[arm][i] > separators[arm]:
                true_PA_sep[arm][i] -= 180


# Поиск пересечения двух множеств.
# Для этого сначала соединим массивы с возможнами
# углами наклона и позиционного угла в один двумерный массив
true_i_PA = []
for i in xrange(nArms):
    true_i_PA.append([])
    for j in xrange(len(true_i[i])):
        true_i_PA[-1].append((true_i[i][j], true_PA[i][j]))

# Найдем пересечение области для всех спиралей
true_i_PA_all = list(reduce(lambda x, y: set(x) & set(y), true_i_PA))
# разобъем его обратно на массив с наклонами и массив с позиционными углами
true_i_all = [i[0] for i in true_i_PA_all]
true_PA_all = [i[1] for i in true_i_PA_all]

true_PA_all_plot = true_PA_all[:] # <- позиционные углы для отображения на рисунке
for i in xrange(len(true_PA_all)):
    if true_PA_all_plot[i] < 0.0:
        true_PA_all_plot[i] += 180

if len(true_i_all) != 0:
    # Если множества имеют пересечение, то вычисляем средний
    # наклон и позиционный угол только по общей части
    incl = mean(true_i_all)
    inclPM = std(true_i_all)
    PA = mean(true_PA_all)
    PAPM = std(true_PA_all)
else:
    incl = 0
    inclPM = 0
    PA = 0
    PAPM = 0
    for i in xrange(nArms):
        incl += mean(true_i[i])
        inclPM += (0.5 * std(true_i[i])) ** 2.0
        PA += mean(true_PA_sep[i])
        PAPM += (0.5 * std(true_PA_sep[i])) ** 2.0
    incl /= nArms
    inclPM **= 0.5
    PA /= nArms
    PAPM **= 0.5

for i in xrange(nArms):
    for j in xrange(len(true_PA_sep[i])):
        if true_PA_sep[i][j] < 0:
            true_PA_sep[i][j] += 180
    
if PA < -0.1:
    PA += 180

if isnan(incl) or isnan(PA):
    print "There is no proper orientation"
    sys.exit(0)

xlabel('$PA  [deg]$', size=14)
ylabel('$i  [deg]$', size=14)
print round(incl, 1), round(PA, 1)

# Отсортируем области по размеру (сначала большие), чтобы
# рисовать сначала самые длинные
srtI = sorted(true_i, key=lambda x: 1./len(x)) 
srtPA = sorted(true_PA_sep, key=lambda x: 1./len(x))

li = [] # лист графиков
color = ['b', 'r', 'g', 'c', 'm', 'y'] # лист цветов
for arm in xrange(nArms):
    li.append(plot(srtPA[arm], srtI[arm], 'ro'))
    setp(li[arm], 'markersize', 5)
    setp(li[arm], 'markerfacecolor', color[arm])
    setp(li[arm], 'markeredgecolor', color[arm])

l = plot(true_PA_all_plot, true_i_all, 'ro')
plot(PA, incl, 'o', color = 'k', markersize=7)
setp(l, 'markersize', 5)
setp(l, 'markerfacecolor', '0.5')
setp(l, 'markeredgecolor', '0.5')
text(130, 80, '$i=%i \degree \pm %i \degree$' % (int(round(incl)), int(round(inclPM))), size=18)
text(130, 75, '$PA=%i \degree \pm %i \degree$' % (int(round(PA)), int(round(PAPM))), size=18)
axis([0.0, 180.0, 0.0, 90.0])
grid(True)
show()
