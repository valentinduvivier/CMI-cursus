# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 16:41:33 2020

@author: vltn0
"""
import matplotlib.pyplot as plt
import numpy as np

stress = np.array([
                  [34215, 1100579.498373958],
                  [27798, 991191.6442907607],
                  [24005, 925641.8971082912],
                  [20845, 872878.2924010159],
                  [18737, 830518.0958818207],
                  [17062, 799429.9687967362],
                  [14507, 747570.6380728119],
                  [12719, 707102.5264564146],
                  [10968, 662390.3453617309],
                  [9551, 628460.0131934556],
                  [7968, 583587.4992575366],
                  [6609, 545837.0390081026],
                  [5678, 513035.7036364893],
                  [5107, 494848.20965236356],
                  [4341, 462288.9796509807],
                  [3713, 440708.36826711806],
                  [3303, 419018.365589503],
                  [2611, 383145.83487306046],
                  [2238, 353580.3088452751],
                  [1969, 345733.84975368227],
                  [1501, 343367.95797897934],
                  [1220, 340612.0595397733],
                  [1030, 338597.9819370444],
                  [866, 330923.792442]
                  ])

stress_2 = np.array([
                  [31289, 1306015.4577835416],
                  [27762, 1301757.1006631001],
                  [22401, 1300289.9672823104],
                  [18611, 1301165.1303069105],
                  [15750, 1297101.6890780944],
                  [11712, 1300933.0472912767],
                  [9686, 1300604.4797397552],
                  [7354, 1299614.7433536802],
                  [5805, 1297537.517257708],
                  [4714, 1296567.9102495268],
                  [3534, 1296360.7724708011],
                  [2902, 1302520.3024280441],
                  [2163, 1298533.419095175],
                  [1672, 1298785.8989027694],
                  [1171, 1292526.4651884406],
                  [876, 1291529.7131127617],
                  [633, 1288039.6159940092],
                  [470, 1289899.5790318872],
                  [324, 1283860.7137985446],
                  [258, 1283972.6663992302],
                  [178, 1279058.8870710498],
                  [135, 1264303.5072904779],
                  [115, 1241617.7772862227],
                  [112, 1274559.5801715425],
                  [108, 1241075.3734157563]
                  ])

# display
fig = plt.figure()

axes1 = fig.add_axes([.1, .1, .8, .8])
#axes2 = fig.add_axes([.3, .4, .5, .3])
axes2 = fig.add_axes([.1, .1, .8, .8])

axes1.plot(stress[:,0], stress[:,1], label = 'varying nods in A')

axes1.set_title('Maximum stress in A over mesh precision')
axes1.set_xlabel('number nods')
axes1.set_ylabel('max stress')

axes1.legend(loc=4)
axes1.grid(True)

# ----------------------------------

axes2.plot(stress_2[:,0], stress_2[:,1], label = 'fix nods in A')

axes2.set_title('Maximum stress in A over mesh precision')
axes2.set_xlabel('number nods')
axes2.set_ylabel('max stress')

axes1.legend(loc=4)
axes1.grid(True)

n = 2000
red_line_x = np.linspace(n, n, np.size(stress_2,0))
red_line_y = np.linspace(0, max(stress_2[:,1])+10**5, np.size(stress_2,0))

plt.plot(red_line_x, red_line_y, '-r', label = f'{n} nods case')

plt.legend(loc = 4)
"""
# insert
axes2.plot(stress[:,0], stress[:,1], 'g-', label = 'T3 element')

axes2.set_title('zommed graph')
axes2.set_xlabel('number nods')
axes2.set_ylabel('max stress')

plt.grid(True)
axes2.legend()

axes2.axis([-150., 1000., 0, 3*10**6])
"""
# ----------------------------------

"""
n = 1000
red_line_x = np.linspace(n, n, np.size(stress,0))
red_line_y = np.linspace(0, max(stress[:,1]), np.size(stress,0))

plt.figure()
plt.plot(stress[:,0], stress[:,1], label = 'T3 element')
plt.plot(red_line_x, red_line_y, '-r', label = '1000 nods case')

plt.title('Maximum stress in A over mesh precision')
plt.xlabel('number nods')
plt.ylabel('max stress')

plt.legend(loc=5)
plt.grid('True')
"""