# -*- coding: utf-8 -*-
"""
Created on Sat Dec  5 19:42:18 2020

@author: vltn0
"""

import numpy as np
import matplotlib.pyplot as plt

disp_T6 = np.array([
                   [28759, 0.39344196558231914],
                   [25576, 0.3932490992052015],
                   [22162, 0.3929456298674559],
                   [17407, 0.3924422463759092],
                   [13372, 0.39177319125877885],
                   [10690, 0.39109451031447384],
                   [7900, 0.3901603285070236],
                   [6616, 0.38956348530639445],
                   [5269, 0.3884758231581368],
                   [4125, 0.3873634874911611],
                   [3291, 0.38611894333105806],
                   [2286, 0.3835917364082317],
                   [1550, 0.38065744884787384],
                   [1008, 0.3756269257323415],
                   [690, 0.36907437965724416],
                   [543, 0.3642260536685872],
                   [240, 0.36772282904648407],
                   [174, 0.36857651026253824],
                   [148, 0.36869955023618667],
                   [130, 0.3681212197791762],
                   [94, 0.35642840640797574],
                   [36, 0.2506673738923565],
                   [16, 0.2511409394190891]
                    ])

n = 1000
red_line_x = np.linspace(n, n, np.size(disp_T6,0))
red_line_y = np.linspace(min(disp_T6[:,1]), max(disp_T6[:,1]), np.size(disp_T6,0))

plt.plot(disp_T6[:,0], disp_T6[:,1])
plt.plot(red_line_x, red_line_y, '-r', label = '1000 nods case')

plt.title('Maximum displacement in C - T6')
plt.xlabel('number nods')
plt.ylabel('max displacement')

plt.grid('True')