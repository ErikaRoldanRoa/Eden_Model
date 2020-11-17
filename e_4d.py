#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 09:55:14 2020

@author: err
"""

import random
import matplotlib.pyplot as plt
from tqdm import tqdm
import gudhi as gd
import numpy as np
import collections
import matplotlib.patches as mpatches
from scipy.optimize import curve_fit

"""GUDHI"""
def convert_gudhi(process):

    min_x = min(process, key=lambda x:x[0])
    max_x = max(process, key=lambda x:x[0])
    min_y = min(process, key=lambda x:x[1])
    max_y = max(process, key=lambda x:x[1])
    min_z = min(process, key=lambda x:x[2])
    max_z = max(process, key=lambda x:x[2])
    min_w = min(process, key=lambda x:x[3])
    max_w = max(process, key=lambda x:x[3])
    dimension = 4
    long = max_x[0] + ((-1)*(min_x[0])) + 1
    wide = max_y[1] + ((-1)*(min_y[1])) + 1
    deep = max_z[2] + ((-1)*(min_z[2])) + 1
    blup = max_w[3] + ((-1)*(min_w[3])) + 1

    time = len(process)
    filename = '4d/gudhi_'+str(time)+'.txt'

    total = long*wide*deep*blup
    pbar = tqdm(total=total)

    with open(filename, 'w+') as f:
        f.writelines('%d\n' % dimension)
        f.writelines('%d\n' % long)
        f.writelines('%d\n' % wide)
        f.writelines('%d\n' % deep)
        f.writelines('%d\n' % blup)
        for w in range(min_w[3], max_w[3] + 1):
            for z in range(min_z[2], max_z[2] + 1):
                for i in range(min_y[1], max_y[1] + 1):
                    for j in range(min_x[0], max_x[0] + 1):
                        pbar.update(1)
                        if (j, i, z, w) in process:
                            f.writelines('%d\n' % process.index((j, i, z, w)))
                        else:
                            f.writelines('inf\n')
    return filename

def gudhi_analysis(filename, final_barcode, time):
    print('What is the minimum length of the interval? Enter 3 numbers (b1, b2, b3) one by one. ')
    length = []
    for i in range(3):
        while True:
            try:
                x = int(input())
                length.append(x)
                break
            except ValueError:
                print("Oops!  That was no valid number.  Try again...")
    print("\nCreating Cubical Complex...")
    eden_model = gd.CubicalComplex(perseus_file=filename)
    eden_model.persistence()
    # A = eden_model.persistence_intervals_in_dimension(1)
    # B = [elem for elem in A if elem[1] == float('inf')]
    barcode_gudhi = eden_model.persistence_intervals_in_dimension(3)
    final = np.array(final_barcode)
    barcode_gudhi_sorted = barcode_gudhi.sort()
    final_sorted = final.sort()
    if barcode_gudhi_sorted == final_sorted:
        print("\nGudhi Barcode agrees with our Barcode!")
    else:
        print("!!!!")

    print("\nDrawing Barcode for Betti_1...")
    pers_1 = [x for x in eden_model.persistence(min_persistence=length[0]) if x[0] == 1]
    fig, ax = plt.subplots()
    gd.plot_persistence_barcode(persistence=pers_1, max_barcodes=1000)
    ax.set_title(r'Persistence Barcode $\beta_1$')
    plt.savefig('4d/'+str(int(time/1000))+'k/barcode_1_'+str(time)+'.png', dpi=1200)

    print("\nDrawing Barcode for Betti_2...")
    pers_2 = [x for x in eden_model.persistence(min_persistence=length[1]) if x[0] == 2]
    fig, ax = plt.subplots()
    gd.plot_persistence_barcode(pers_2, max_barcodes=1000)
    ax.set_title(r'Persistence Barcode $\beta_2$')
    plt.savefig('4d/'+str(int(time/1000))+'k/barcode_2_'+str(time)+'.png', dpi=1200)

    print("\nDrawing Barcode for Betti_3...")
    pers_3 = [x for x in eden_model.persistence(min_persistence=length[2]) if x[0] == 3]
    fig, ax = plt.subplots()
    gd.plot_persistence_barcode(pers_3, max_barcodes=1000)
    ax.set_title(r'Persistence Barcode $\beta_3$')
    plt.savefig('4d/'+str(int(time/1000))+'k/barcode_3_'+str(time)+'.png', dpi=1200)

def convert_perseus_2(Process):

    dimension = 4
    with open('1000000_1_4D.txt','w') as f:
        f.writelines( '%d\n' % dimension)
        i = 0
        for x in Process:
            i = i + 1
            y = (x[0], x[1], x [2], x[3], i)
            f.writelines( '%s %s %s %s %s\n' % y)

"""GROWING"""
def grow_eden(t):

    process = [(0, 0, 0, 0)]
    perimeter_len = [8]
    eden, perimeter = start_eden()

    l = len(perimeter)

    holes = {}
    total_holes = 0
    barcode = {}
    tags = []
    created_holes = []

    betti_3_total = 0
    betti_3_vector = [0]
    betti_3_total_vector = [0]

    pbar = tqdm(total=t)
    pbar.update(1)
    for i in range(1, t):
        pbar.update(1)
        l = len(perimeter)
        x = random.randint(0, l-1)
        tile_selected = perimeter[x]
        perimeter.pop(x)

        eden[tile_selected][0] = 1
        process = process + [tile_selected]

        eden, perimeter, nearest_n, nearest_n_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        betti_3, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_3(eden, tile_selected, nearest_n, nearest_n_tiles, barcode, i, holes, total_holes, created_holes, tags)
        betti_3_vector = betti_3_vector + [betti_3]
        betti_3_total = betti_3_total + betti_3
        betti_3_total_vector += [betti_3_total]

        l = len(perimeter)
        perimeter_len = perimeter_len + [l]
    final_barcode = barcode_forest(barcode, tags)

    return eden, perimeter, betti_3_vector, barcode, holes, betti_3_total, created_holes, process,\
           perimeter_len, betti_3_total_vector, final_barcode

def grow_eden_debugging(t, ordered_tiles):

    eden, perimeter = start_eden()
    perimeter_len = [8]

    holes = {}
    total_holes = 0
    barcode = {}
    tags = []
    created_holes = []

    betti_3_total = 0
    betti_3_vector = [0]
    betti_3_total_vector = []

    pbar = tqdm(total=t)
    pbar.update(1)

    for i in (range(1, t)):
        pbar.update(1)

        tile_selected = ordered_tiles[i][:4]
        perimeter.remove(tile_selected)
        eden[tile_selected][0] = 1

        l = len(perimeter)
        perimeter_len = perimeter_len + [l]

        eden, perimeter, nearest_n, nearest_n_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        perimeter_len = perimeter_len + [len(perimeter)]

        betti_3, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_3(eden, tile_selected, nearest_n, nearest_n_tiles, barcode, i, holes, total_holes, created_holes, tags)
        betti_3_vector = betti_3_vector + [betti_3]
        betti_3_total = betti_3_total + betti_3
        betti_3_total_vector += [betti_3_total]
    final_barcode = barcode_forest(barcode, tags)

    return eden, perimeter, betti_3_vector, barcode, holes, betti_3_total, created_holes,\
           perimeter_len, betti_3_total_vector, final_barcode

"""PLOTTING"""
def draw_frequencies_3(dict, time, changes):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    l = len(dict[0])

    ch_5 = [i for i, j in enumerate(changes) if j == 5]
    y_5 = []
    for x in ch_5:
        y_5 += [dict[5][x+1]]

    sh = []
    for j in np.arange(-1, 3):
        sh.append(next((i for i, x in enumerate(dict[j]) if x), 0))
    shift = max(sh)

    ax.plot(range(shift, l), dict[-1][shift:], color='tab:olive', label='-1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[0][shift:], color='black', label='0',  linewidth=0.75)
    ax.plot(range(shift, l), dict[1][shift:], color='tab:red', label='+1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[2][shift:], color='tab:orange', label='+2',  linewidth=0.75)
    ax.plot(range(shift, l), dict[3][shift:], color='tab:green', label='+3',  linewidth=0.75)
    ax.plot(range(shift, l), dict[4][shift:], color='tab:blue', label='+4',  linewidth=0.75)
    # shift = next((i for i, x in enumerate(dict[3]) if x), 0)
    # ax.plot(range(shift, l), dict[3][shift:], color='tab:purple', label='+3',  linewidth=0.75)
    # ax.plot(range(shift, l), dict[4][shift:], color='tab:brown', label='4',  linewidth=0.75)
    plt.scatter(ch_5, y_5, s=5, marker='o', color="tab:brown", label='+5')
    # print(ch_5)
    # plt.scatter(ch_m_4, y_m_4, s = 10, marker='o', color="black", label='-4')

    plt.yscale('log')
    # ax.set_title('betti_1 frequencies')
    ax.set_ylabel(r'Frequency of Change in $\beta_3$')
    ax.set_xlabel('t')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # ax.ticklabel_format(useOffset=False)
    ax.legend(loc=1, prop={'size': 6})
    fig.savefig('4d/'+str(int(time/1000))+'k/fr_b_3_'+str(time)+'.png', format='png', dpi=1200)
    # plt.show()
    plt.close()

def plot_b_per(Betti_3_total_vector, Per, time, N):
    n = int(time/10)

    def func(x, a, b):
        return a * x ** b
    ydata_f = Betti_3_total_vector
    xdata_f = range(len(ydata_f))
    ydata = ydata_f[N:]
    xdata = xdata_f[N:]
    plt.xscale('log')
    plt.yscale('log')
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(xdata_f[n:], ydata_f[n:], 'm-', label=r'$\beta_3(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)

    plt.plot(xdata_f[n:], func(xdata_f[n:], *popt), 'm--', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt), linewidth=0.75)

    ydata = Per
    xdata = range(len(ydata))
    plt.plot(xdata[n:], ydata[n:], color='orange', linestyle='solid', label=r'$P(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata[n:], func(xdata[n:], *popt), color='orange', linestyle='dashed', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt),  linewidth=0.75)

    plt.xlabel('t')
    plt.ylabel('data')
    plt.legend(prop={'size': 6}, loc=2)
    plt.tight_layout()

    plt.savefig('4d/'+str(int(time/1000))+'k/per-b-time_'+str(time)+'.png', dpi=1200)
    # plt.show()
    plt.close()


"""SUPPLEMENTARY FUNCTIONS"""
def read_eden_txt(filename):
    eden = []
    for t in open(filename).read().split('), ('):
        # print(t)
        a, b, c, d, t = t.strip('()[]').split(',')
        a = a.strip()
        b = b.strip()
        c = c.strip()
        d = d.strip(')')
        t = t.strip(")]\n")
        eden.append(((int(a), int(b), int(c), int(d)), float(t)))
    return eden

def num_holes(created_holes, holes):
    tricube = dict.fromkeys(['l','i'], 0)
    tricube_f = dict.fromkeys(['l','i'], 0)
    tetracube = dict.fromkeys(['I', 'L', 'T', 'O', 'Z', 'A1', 'A2'], 0)
    tetracube_f = dict.fromkeys(['I', 'L', 'T', 'O', 'Z', 'A1', 'A2'], 0)

    total_3 = [i for i in created_holes if i[-1] == 3]
    final_3 = [holes[i] for i in holes if len(holes[i]) == 3]
    total_4 = [i for i in created_holes if i[-1] == 4]
    final_4 = [holes[i] for i in holes if len(holes[i]) == 4]

    # tricube case
    for x in total_3:
        j = 0
        for i in range(3):
            if x[-2][0][i] == x[-2][1][i] == x[-2][2][i]:
                j += 1
        long = j == 2
        tricube['i'] += long
        tricube['l'] += (1 - long)
    for x in final_3:
        j = 0
        for i in range(3):
            if x[0][i] == x[1][i] == x[2][i]:
                j += 1
        long = j == 2
        tricube_f['i'] += long
        tricube_f['l'] += (1 - long)

    # tetracube case
    for x in total_4:
        dist = distances(x[-2])
        if dist == [2, 2, 3]:
            typ = 'I'
        elif dist == [1.4, 2, 2.2]:
            typ = 'L'
        elif dist == [1.4, 1.4, 2]:
            typ = 'T'
        elif dist == [1, 1.4, 1.4]:
            typ = 'O'
        elif dist == [1.4, 1.4, 2.2]:
            typ = 'Z'
        elif dist == [1.4, 1.4, 1.4]:
            typ = 'A2'
        else:
            typ = 'A1'
        tetracube[typ] += 1
    for x in final_4:
        dist = distances(x)
        if dist == [2,2,3]:
            typ = 'I'
        elif dist == [1.4, 2, 2.2]:
            typ = 'L'
        elif dist == [1.4, 1.4, 2]:
            typ = 'T'
        elif dist == [1, 1.4, 1.4]:
            typ = 'O'
        elif dist == [1.4, 1.4, 2.2]:
            typ = 'Z'
        elif dist == [1.4, 1.4, 1.4]:
            typ = 'A2'
        else:
            typ = 'A1'
        tetracube_f[typ] += 1
    return tricube, tricube_f, tetracube, tetracube_f

def distances(hole):
    hole = [np.array(k) for k in hole]
    dist = []
    for i in range(4):
        for j in range(4):
            if i < j:
                dist.append(np.linalg.norm(hole[i]-hole[j]))
    dist.sort()
    dist = [round(i, 1) for i in dist]
    return dist[3:]

def return_frequencies_3(vect, time):
    changes = [vect[i+1]-vect[i] for i in range(len(vect)-1)]
    # values = list(set(changes))
    # values.sort()
    values = [-1, 0, 1, 2, 3, 4, 5]
    freq = {i: [0] for i in values}

    for i in tqdm(range(1, time+1), position=0, leave=True):
        counter = collections.Counter(changes[:i])
        for k in values:
            freq[k].append(counter[k]/i)
    return freq, changes

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def Diff(li1, li2):
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

def start_eden():
    Eden = {(0, 0, 0, 0) : [1, 0, 0],
            (0, 0, 1, 0) : [0, 1, 0],
            (0, 0, -1, 0): [0, 1, 0],
            (1, 0, 0, 0) : [0, 1, 0],
            (-1, 0, 0, 0): [0, 1, 0],
            (0, 1, 0, 0) : [0, 1, 0],
            (0,-1, 0, 0) : [0, 1, 0],
            (0, 0, 0, 1) : [0, 1, 0],
            (0, 0, 0, -1): [0, 1, 0]}
    Perimeter = [(0, 0, 1, 0), (0, 0, -1, 0), (1, 0, 0, 0), (-1, 0, 0, 0), (0, 1, 0, 0), (0, -1, 0, 0), (0, 0, 0, 1), (0, 0, 0, -1)]
    return Eden, Perimeter

def actualize_neighbors(tile_selected, Eden, Perimeter):

    n3 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2], tile_selected[3]]
    n4 = [tile_selected[0]  -1, tile_selected[1], tile_selected[2], tile_selected[3]]
    n5 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2], tile_selected[3]]
    n6 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2], tile_selected[3]]
    n1 = [tile_selected[0], tile_selected[1], tile_selected[2] + 1, tile_selected[3]]
    n2 = [tile_selected[0], tile_selected[1], tile_selected[2] - 1, tile_selected[3]]
    n7 = [tile_selected[0], tile_selected[1], tile_selected[2], tile_selected[3] + 1]
    n8 = [tile_selected[0], tile_selected[1], tile_selected[2], tile_selected[3] - 1]

    n1 = tuple(n1)
    n2 = tuple(n2)
    n3 = tuple(n3)
    n4 = tuple(n4)
    n5 = tuple(n5)
    n6 = tuple(n6)
    n7 = tuple(n7)
    n8 = tuple(n8)

    nearest_n_tiles = [n1, n2, n3, n4, n5, n6, n7, n8]

    nearest_n = [0, 0, 0, 0, 0, 0, 0, 0]

    if n1 in Eden:
        Eden[n1][1] = Eden[n1][1] + 1
        if Eden[n1][0] == 1:
            nearest_n[0] = 1
    else:
        Eden[n1] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n1]

    if n2 in Eden:
        Eden[n2][1] = Eden[n2][1] + 1
        if Eden[n2][0] == 1:
            nearest_n[1] = 1
    else:
        Eden[n2] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n2]

    if n3 in Eden:
        Eden[n3][1] = Eden[n3][1] + 1
        if Eden[n3][0] == 1:
            nearest_n[2] = 1
    else:
        Eden[n3] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n3]

    if n4 in Eden:
        Eden[n4][1] = Eden[n4][1] + 1
        if Eden[n4][0] == 1:
            nearest_n[3] = 1
    else:
        Eden[n4] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n4]

    if n5 in Eden:
        Eden[n5][1] = Eden[n5][1] + 1
        if Eden[n5][0] == 1:
            nearest_n[4] = 1
    else:
        Eden[n5] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n5]

    if n6 in Eden:
        Eden[n6][1] = Eden[n6][1] + 1
        if Eden[n6][0] == 1:
            nearest_n[5] = 1
    else:
        Eden[n6] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n6]

    if n7 in Eden:
        Eden[n7][1] = Eden[n7][1] + 1
        if Eden[n7][0] == 1:
            nearest_n[6] = 1
    else:
        Eden[n7] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n7]

    if n8 in Eden:
        Eden[n8][1] = Eden[n8][1] + 1
        if Eden[n8][0] == 1:
            nearest_n[7] = 1
    else:
        Eden[n8] = [0, 1, Eden[tile_selected][2]]
        Perimeter = Perimeter + [n8]

    return Eden, Perimeter, nearest_n, nearest_n_tiles

def add_neighbours_bds(bds, j, iterations, num_possible_components, merged, finished, Eden):

    tile_selected = bds[j][iterations]


    n3 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2], tile_selected[3]]
    n4 = [tile_selected[0]  -1, tile_selected[1], tile_selected[2], tile_selected[3]]
    n5 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2], tile_selected[3]]
    n6 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2], tile_selected[3]]
    n1 = [tile_selected[0], tile_selected[1], tile_selected[2] + 1, tile_selected[3]]
    n2 = [tile_selected[0], tile_selected[1], tile_selected[2] - 1, tile_selected[3]]
    n7 = [tile_selected[0], tile_selected[1], tile_selected[2], tile_selected[3] + 1]
    n8 = [tile_selected[0], tile_selected[1], tile_selected[2], tile_selected[3] - 1]


    n1 = tuple(n1)
    n2 = tuple(n2)
    n3 = tuple(n3)
    n4 = tuple(n4)
    n5 = tuple(n5)
    n6 = tuple(n6)
    n7 = tuple(n7)
    n8 = tuple(n8)

    nearest_n_tiles = [n1, n2, n3, n4, n5, n6, n7, n8]

    nearest_n = [0, 0, 0, 0, 0, 0, 0, 0]

    if n1 in Eden:
        if Eden[n1][0] == 0:
            nearest_n[0] = 1
    else:
        nearest_n[0] = 1

    if n2 in Eden:
        if Eden[n2][0] == 0:
            nearest_n[1] = 1
    else:
        nearest_n[1] = 1

    if n3 in Eden:
        if Eden[n3][0] == 0:
            nearest_n[2] = 1
    else:
        nearest_n[2] = 1

    if n4 in Eden:
        if Eden[n4][0] == 0:
            nearest_n[3] = 1
    else:
        nearest_n[3] = 1

    if n5 in Eden:
        if Eden[n5][0] == 0:
            nearest_n[4] = 1
    else:
        nearest_n[4] = 1

    if n6 in Eden:
        if Eden[n6][0] == 0:
            nearest_n[5] = 1
    else:
        nearest_n[5] = 1

    if n7 in Eden:
        if Eden[n7][0] == 0:
            nearest_n[6] = 1
    else:
        nearest_n[6] = 1

    if n8 in Eden:
        if Eden[n8][0] == 0:
            nearest_n[7] = 1
    else:
        nearest_n[7] = 1

    for i in range(0,8):
        if nearest_n[i] == 1:
            if nearest_n_tiles[i] not in bds[j]:
                bds[j] = bds[j] + [nearest_n_tiles[i]]
            for t in range(0, num_possible_components):
                if nearest_n_tiles[i] in bds[t]:
                    if t < j:
                        merged[j] = 1
                        finished[j] = 1
                    if t > j:
                        merged[t] = 1
                        finished[t] = 1
    return bds, merged, finished

def increment_betti_3(Eden, tile_selected, nearest_n, nearest_n_tiles, barcode, time, Holes, total_holes, Created_Holes, tags):
    if Eden[tile_selected][2] == 0:
        per = 1                       # This is 1 if the tile added was in the out perimeter
    else:
        num_hole = Eden[tile_selected][2]
        per = 0
  # In this case the tile added was in a hole

    betti_3 = 0

    if sum(nearest_n) == 8:
        betti_3 = - 1
        barcode[num_hole][1][1] = float(time + 2)   #Are we covering a hole that was never divided?
        Holes[num_hole].remove(tile_selected)
        if Holes[num_hole] == []:
            Holes.pop(num_hole)

    if sum(nearest_n) == 7:
        betti_3 = 0
        if per == 0:
            Holes[num_hole].remove(tile_selected)
    #print(nearest_n)
    if sum(nearest_n) != 8 and sum(nearest_n) != 7:
        num_posible_components = 0
        bds = []
        iterations = 0
        for i in range(0,8):
            if nearest_n[i] == 0:
                num_posible_components = num_posible_components + 1
                bds = bds + [[nearest_n_tiles[i]]]

        finished = [0] * num_posible_components
        merged = finished.copy()


        while sum(finished) < num_posible_components - per:
            for j in range(0, num_posible_components):
                if finished[j] == 0:
                    bds, merged, finished = add_neighbours_bds(bds, j, iterations, num_posible_components, merged, finished, Eden)
                    if (iterations + 1) == len(bds[j]):
                        finished[j] = 1
            iterations = iterations + 1

        betti_3 = (num_posible_components - 1) - sum(merged)
        #print(betti_2, per)
        #At this point we have the bds components and the ones that were not merged will become the holes.
        #Here we actualize the holes and we actualize Hole No in Eden.
        if betti_3 == 0:
            if per == 0:
                Holes[num_hole].remove(tile_selected)

        else:
            if per == 1:
                for i in range(0, num_posible_components):
                    if finished[i] == 1 and merged[i] == 0:
                        total_holes = total_holes + 1
                        Holes[total_holes] = bds[i].copy()

                        for x in bds[i]:
                            if x in Eden:
                                Eden[x][2] = total_holes

                        barcode[total_holes] = [0, [float(time + 2), float(0)], [total_holes]]
                        Created_Holes = Created_Holes + [[barcode[total_holes][2], bds[i].copy(), len(bds[i])]]

            else:
                if barcode[num_hole][0] == 0:
                    tags = tags + [num_hole]
                    barcode[num_hole][0] = 1

                Holes.pop(num_hole)

                for i in range(0, num_posible_components):
                    if  finished[i] == 1 and merged[i] == 0:
                        total_holes = total_holes + 1
                        Holes[total_holes] = bds[i].copy()
                        for x in bds[i]:
                            if x in Eden:
                                Eden[x][2] = total_holes
                        barcode[total_holes] = [1, [float(time + 2), float(0)], barcode[num_hole][2] + [total_holes]]
                        Created_Holes = Created_Holes + [[barcode[total_holes][2], bds[i].copy(), len(bds[i])]]


    return  betti_3, total_holes, Eden, barcode, Holes, Created_Holes, tags

def barcode_forest(barcode, tags):
    bars_pure = []
    bars_hole = []
    for x in barcode:
        if barcode[x][0] == 0:
            bars_pure = bars_pure + [barcode[x][1]]

    for x in tags:
        b = {}
        for elem in barcode:
            if barcode[elem][2][0] == x:
                b[tuple(barcode[elem][2])] = barcode[elem][1]
#        print(b)
        bars_hole = bars_hole + bars_from_tree(b, x)
    return bars_pure + bars_hole

def bars_from_tree(b, tag):
    n = max(len(x) for x in b)
    bars = []
    while n > 0:
        leaves_parent = [ x for x in b if len(x)== n - 1 ]

#        print(leaves_parent)
        possible_leaves = [ x for x in b if len(x)== n ]
#        print(possible_leaves)
        for j in leaves_parent:
            leaves = []
            for x in possible_leaves:
#            print(x)
                root = list(x)
                del root[-1]
                root = tuple(root)
                if  hamming2(j, root) == 0:
                    leaves = leaves + [x]
#            Diff(possible_leaves, leaves)
#            print(leaves)
            if len(leaves) > 0 :
                times = []
                for x in leaves:
                    times = times  + [b[x][1]]
                if 0 in times:
                    ind = times.index(0)
                    for i in range(0, len(leaves)):
                        if i != ind:
                            bars = bars + [b[leaves[i]]]
                else:
                    ind = times.index(max(times))
                    for i in range(0, len(leaves)):
                        if i != ind:
                            bars = bars + [b[leaves[i]] ]
                    b[j][1] = max(times)
        n = n - 1
    bars = bars + [b[(tag,)]]
    return bars




