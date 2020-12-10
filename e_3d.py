#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 21 19:34:33 2019

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
def convert_gudhi(process, folder_name):
    min_x = min(process, key=lambda x: x[0])
    max_x = max(process, key=lambda x: x[0])
    min_y = min(process, key=lambda x: x[1])
    max_y = max(process, key=lambda x: x[1])
    min_z = min(process, key=lambda x: x[2])
    max_z = max(process, key=lambda x: x[2])
    dimension = 3
    long = max_x[0] + ((-1) * (min_x[0])) + 1
    wide = max_y[1] + ((-1) * (min_y[1])) + 1
    deep = max_z[2] + ((-1) * (min_z[2])) + 1
    time = len(process)
    filename = folder_name+'/gudhi_'+str(time)+'.txt'

    total = long*wide*deep
    pbar = tqdm(total=total)

    with open(filename, 'w+') as f:
        f.writelines('%d\n' % dimension)
        f.writelines('%d\n' % long)
        f.writelines('%d\n' % wide)
        f.writelines('%d\n' % deep)
        for z in range(min_z[2], max_z[2] + 1):
            for i in range(min_y[1], max_y[1] + 1):
                for j in range(min_x[0], max_x[0] + 1):
                    pbar.update(1)
                    if (j, i, z) in process:
                        f.writelines('%d\n' % process.index((j, i, z)))
                    else:
                        f.writelines('inf\n')
    pbar.close()
    return filename

def gudhi_analysis(filename, final_barcode, time):
    print('What is the minimum length of the interval? Enter 3 numbers one by one. ')
    length = []
    for i in range(2):
        print("Minimal length for Betti_"+str(i+1)+':')
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
    A = eden_model.persistence_intervals_in_dimension(2)
    B = [elem for elem in A if elem[1] == float('inf')]
    final = np.array(final_barcode)
    A_sorted = A.sort()
    final_sorted = final.sort()
    if A_sorted == final_sorted:
        print("Gudhi Barcode agrees with our Barcode!")
    pers_1 = [x for x in eden_model.persistence(min_persistence=length[0]) if x[0] == 1]
    fig, ax = plt.subplots()
    gd.plot_persistence_barcode(persistence=pers_1, max_barcodes=1000)
    ax.set_title(r'Persistence Barcode $\beta_1$')
    plt.savefig(folder_name+'/barcode_1.png', dpi=1200)

    pers_2 = [x for x in eden_model.persistence(min_persistence=length[1]) if x[0] == 2]
    fig, ax = plt.subplots()
    gd.plot_persistence_barcode(pers_2, max_barcodes=1000)
    ax.set_title(r'Persistence Barcode $\beta_2$')
    plt.savefig(folder_name+'/barcode_2.png', dpi=1200)

def convert_perseus_2(process):
    dimension = 3
    with open('1000000_3D_12_final.txt', 'w') as f:
        f.writelines('%d\n' % dimension)
        i = 0
        for x in process:
            i = i + 1
            y = (x[0], x[1], x[2], i)
            f.writelines('%s %s %s %s\n' % y)

"""GROWING"""
def grow_eden(t):
    vertices = 8
    edges = 12
    faces = 6

    process = [(0, 0, 0)]
    perimeter_len = [6]
    eden, perimeter = start_eden()

    holes = {}
    total_holes = 0
    barcode = {}
    tags = []
    created_holes = []

    betti_2_total = 0
    betti_2_vector_changes = [0]
    betti_2_total_vector = [0]

    betti_1_total_vector = [0]

    skipped = 0
    size = 1

    pbar = tqdm(total=t)
    pbar.update(1)
    # for i in tqdm(range(1, t)):
    # for i in (range(1, t)):
    while size < t:
        l = len(perimeter)
        x = random.randint(0, l - 1)
        tile_selected = perimeter[x]
        perimeter.pop(x)

        eden[tile_selected][0] = 1
        process = process + [tile_selected]

        eden, perimeter, nearest_n, nearest_n_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        n = neighbours(eden, tile_selected)

        betti_2, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_2(eden, tile_selected,
                                                                                            nearest_n, nearest_n_tiles,
                                                                                            barcode, size, holes,
                                                                                            total_holes, created_holes,
                                                                                            tags)

        vertices, edges, faces, v_new, e_new, f_new = actualize_vef(vertices, edges, faces, nearest_n, n)

        pbar.update(1)
        euler_character = euler_characteristic(vertices, edges, faces, size + 1)
        size += 1
        betti_2_vector_changes += [betti_2]
        betti_2_total += betti_2
        betti_2_total_vector += [betti_2_total]

        betti_1_total = return_betti_1(betti_2_total, euler_character)
        betti_1_total_vector += [betti_1_total]

        l = len(perimeter)
        perimeter_len = perimeter_len + [l]

    final_barcode = barcode_forest(barcode, tags)

    pbar.close()
    return eden, perimeter, betti_2_total_vector, betti_2_vector_changes, barcode, holes, betti_1_total, betti_1_total_vector, \
           created_holes, process, perimeter_len, skipped, size,\
           final_barcode

def grow_eden_debugging(t, ordered_tiles):
    vertices = 8
    edges = 12
    faces = 6
    #
    # process = [(0, 0, 0)]

    eden, perimeter = start_eden()

    holes = {}
    total_holes = 0
    barcode = {}
    tags = []
    created_holes = []
    perimeter_len = [6]

    betti_2_total = 0
    betti_2_vector_changes = [0]
    betti_2_total_vector = [0]

    betti_1_total_vector = [0]
    pbar = tqdm(total=t, position=0, leave=True)
    pbar.update(1)
    for i in (range(1, t)):
        pbar.update(1)

        tile_selected = ordered_tiles[i]
        perimeter.remove(tile_selected)
        eden[tile_selected][0] = 1

        eden, perimeter, nearest_n, nearest_n_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        n = neighbours(eden, tile_selected)

        vertices, edges, faces, v, e, f = actualize_vef(vertices, edges, faces, nearest_n, n)
        betti_2, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_2(eden,
                                                                                            tile_selected, nearest_n,
                                                                                            nearest_n_tiles,
                                                                                            barcode, i, holes,
                                                                                            total_holes, created_holes,
                                                                                            tags)

        betti_2_total = betti_2_total + betti_2
        betti_2_total_vector += [betti_2_total]

        euler_character = euler_characteristic(vertices, edges, faces, i)
        betti_1_total = return_betti_1(betti_2_total, euler_character)
        betti_1_total_vector += [betti_1_total]

        l = len(perimeter)
        perimeter_len = perimeter_len + [l]
    pbar.close()
    final_barcode = barcode_forest(barcode, tags)

    return eden, perimeter, betti_2_total_vector, betti_1_total_vector, barcode, holes, betti_2_total, betti_1_total, \
           created_holes, perimeter_len, final_barcode  # , tags, final_barcode

"""PLOTTING"""
def draw_frequencies_1(dict, changes, folder_name):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    l = len(dict[0])

    ch_4 = [i for i, j in enumerate(changes) if j == 4]
    y_4 = []
    for x in ch_4:
        y_4 += [dict[4][x+1]]

    sh = []
    for j in np.arange(-2, 3):
        sh.append(next((i for i, x in enumerate(dict[j]) if x), 0))
    shift = max(sh)

    shiftt = next((i for i, x in enumerate(dict[-3]) if x), 0)
    ax.plot(range(shiftt, l), dict[-3][shiftt:], color='tab:olive', label='-3',  linewidth=0.75)
    ax.plot(range(shift, l), dict[-2][shift:], color='black', label='-2',  linewidth=0.75)
    ax.plot(range(shift, l), dict[-1][shift:], color='tab:red', label='-1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[0][shift:], color='tab:orange', label='0',  linewidth=0.75)
    ax.plot(range(shift, l), dict[1][shift:], color='tab:green', label='+1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[2][shift:], color='tab:blue', label='+2',  linewidth=0.75)
    shift = next((i for i, x in enumerate(dict[3]) if x), 0)
    ax.plot(range(shift, l), dict[3][shift:], color='tab:purple', label='+3',  linewidth=0.75)
    # ax.plot(range(shift, l), dict[4][shift:], color='tab:brown', label='4',  linewidth=0.75)
    plt.scatter(ch_4, y_4, s=5, marker='o', color="tab:brown", label='+4')
    # plt.scatter(ch_m_4, y_m_4, s = 10, marker='o', color="black", label='-4')

    plt.yscale('log')
    # ax.set_title('betti_1 frequencies')
    ax.set_ylabel(r'Frequency of Change in $\beta_1$')
    ax.set_xlabel('t')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # ax.ticklabel_format(useOffset=False)
    ax.legend(loc=1, prop={'size': 6})
    fig.savefig(folder_name+'/fr_b_1.png', format='png', dpi=1200)
    # plt.show()
    plt.close()

def draw_frequencies_2(dict, changes, folder_name):
    # fig = plt.figure()
    fig, ax = plt.subplots()
    l = len(dict[0])

    ch_4 = [i for i, j in enumerate(changes) if j == 4]
    y_4 = []
    for x in ch_4:
        y_4 += [dict[4][x+1]]

    sh = []
    for j in np.arange(-1, 3):
        sh.append(next((i for i, x in enumerate(dict[j]) if x), 0))
    shift = max(sh)

    ax.plot(range(shift, l), dict[-1][shift:], color='tab:red', label='-1', linewidth=0.75)
    ax.plot(range(shift, l), dict[0][shift:], color='tab:orange', label='0', linewidth=0.75)
    ax.plot(range(shift, l), dict[1][shift:], color='tab:green', label='+1', linewidth=0.75)
    ax.plot(range(shift, l), dict[2][shift:], color='tab:blue', label='+2', linewidth=0.75)
    shift = next((i for i, x in enumerate(dict[3]) if x), 0)
    ax.plot(range(shift, l), dict[3][shift:], color='tab:purple', label='+3', linewidth=0.75)
    ax.scatter(ch_4, y_4, s=5, marker='o', color="tab:brown", label='+4')

    # ax.plot(range(shift, l), dict[4][shift:], color='tab:brown', label='4', linewidth=0.75)

    plt.yscale('log')
    ax.set_ylabel(r'Frequency of Change in $\beta_2$')
    ax.set_xlabel('t')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # ax.ticklabel_format(useOffset=False)
    plt.legend(loc=1, prop={'size': 6})
    fig.savefig(folder_name+'/fr_b_2.png', format='png', dpi=1200)
    # plt.show()
    plt.close()

def draw_tri_tetra(tri, tri_f, tetra, tetra_f, folder_name):
    width = 0.35
    labels = list(tri)+list(tetra)
    x = np.arange(len(labels))
    fig, ax = plt.subplots()
    plt.yscale('log')

    er = [0,0,0,0]
    try:
        ax.bar(x[:2]-width/2, list(tri.values()), width, label='Trominoes Total', color='navy')
    except ValueError:
        er[0] = 1
    try:
        ax.bar(x[2:]-width/2, list(tetra.values()), width, label='Tetrominoes Total', color='chocolate')
    except ValueError:
        er[2] = 1
    try:
        ax.bar(x[:2]+width/2, tri_f.values(), width, label='Trominoes Final', color='royalblue')
    except ValueError:
        er[1] = 1
    try:
        ax.bar(x[2:]+width/2, tetra_f.values(), width, label='Tetrominoes Final', color='orange')

    except ValueError:
        er[3] = 1
        # ax.set_xticks(x)
        # ax.set_xticklabels(labels)
        # handles, labels = ax.get_legend_handles_labels()
        # # labels = [labels[0], labels[2], labels[1]]
        # patch = mpatches.Patch(color='orange', label='Tetrominoes Final', linewidth=0.35)
        # handles.append(patch)
        # labels.append('Tetrominoes Final')
        # ax.legend(handles, labels, loc='upper right')
    if er == [0, 0, 0, 0]:
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper right')
    else:
        if er == [1,1,1,1]:
            handles = []
        else:
            handles, labelss = ax.get_legend_handles_labels()
        if er[0] == 1:
            # handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color='navy', label='Trominoes Total', linewidth=0.35)
            handles.append(patch)
        if er[1] == 1:
            # handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color='royalblue', label='Trominoes Final', linewidth=0.35)
            handles.append(patch)
        if er[2] == 1:
            # handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color='chocolate', label='Tetrominoes Total', linewidth=0.35)
            handles.append(patch)
        if er[3] == 1:
            # handles, labels = ax.get_legend_handles_labels()
            patch = mpatches.Patch(color='orange', label='Tetrominoes Final', linewidth=0.35)
            handles.append(patch)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(handles=handles, loc='upper right')

    # Add some text for labels, title and custom x-axis tick labels, etc.
    ax.set_ylabel('Frequency of Number of Holes')
    ax.set_xlabel('Type of a Hole')
    fig.tight_layout()
    fig.savefig(folder_name+'/tri-tetra-cubes.png', format='png', dpi=1200)
    # plt.show()
    plt.close()

def plot_b_per(Betti_1_total_vector, Betti_2_total_vector, Per, time, N, folder_name):
    n = int(time/10)

    def func(x, a, b):
        return a * x ** b
    ydata_f = Betti_1_total_vector
    xdata_f = range(len(ydata_f))
    ydata = ydata_f[N:]
    xdata = xdata_f[N:]
    plt.xscale('log')
    plt.yscale('log')
    # plt.ticklabel_format(style='sci', axis='x', scilimits=(0,0))
    # plt.ticklabel_format(style='sci', axis='y', scilimits=(0,0))
    plt.plot(xdata_f[n:], ydata_f[n:], 'm-', label=r'$\beta_1(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)

    plt.plot(xdata_f[n:], func(xdata_f[n:], *popt), 'm--', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt), linewidth=0.75)

    ydata_f = Betti_2_total_vector
    xdata_f = range(len(ydata_f))
    ydata = ydata_f[N:]
    xdata = xdata_f[N:]
    plt.plot(xdata_f[n:], ydata_f[n:], 'b-', label=r'$\beta_2(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata_f[n:], func(xdata_f[n:], *popt), 'b--', label=r'fit: $y=%5.3f x^{%5.3f}$' % tuple(popt),  linewidth=0.75)

    ydata = Per
    xdata = range(len(ydata))
    plt.plot(xdata[n:], ydata[n:], color='orange', linestyle='solid', label=r'$P(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata[n:], func(xdata[n:], *popt), color='orange', linestyle='dashed', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt),  linewidth=0.75)

    plt.xlabel('t')
    plt.ylabel('data')
    plt.legend(prop={'size': 6}, loc=2)
    plt.tight_layout()

    plt.savefig(folder_name+'/per-b-time.png', dpi=1200)
    # plt.show()
    plt.close()


"""SUPPLEMENTARY FUNCTIONS"""
def read_eden_txt(filename):
    eden = []
    for t in open(filename).read().split('), ('):
        # print(t)
        a, b, c, t = t.strip('()[]').split(',')
        a = a.strip()
        b = b.strip()
        c = c.strip(')')
        t = t.strip(")]\n")
        eden.append(((int(a), int(b), int(c)), float(t)))
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

def return_frequencies_1(vect, time):
    changes = [vect[i+1]-vect[i] for i in range(len(vect)-1)]
    # values = list(set(changes))
    # values.sort()
    values = [-3, -2, -1, 0, 1, 2, 3, 4]
    freq = {i: [0] for i in values}

    for i in tqdm(range(1, time+1), position=0, leave=True):
        counter = collections.Counter(changes[:i])
        for k in values:
            freq[k].append(counter[k]/i)
    return freq, changes

def return_frequencies_2(vect, time):
    changes = [vect[i+1]-vect[i] for i in range(len(vect)-1)]
    # values = list(set(changes))
    # values.sort()
    values = [-1, 0, 1, 2, 3, 4]
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

def draw_barcode(barcode, time):
    """ """
    fig = plt.figure()
    plt.style.use('ggplot')
    # plt.axis('off')
    # plt.rc('grid', linestyle="-", color='black')
    plt.grid(True)
    plt.rc('grid', linestyle="-", color='gray')
    plt.yticks([])
    plt.gca().set_aspect('equal', adjustable='box')
    i = 0
    for x in barcode:
        if barcode[x][1] == 0:
            plt.plot([barcode[x][0], time], [i, i], 'k-', lw=2)
        else:
            plt.plot([barcode[x][0], barcode[x][1]], [i, i], 'k-', lw=2)
        i = i + 40
    fig.savefig('5000.png')
    # plt.show()

def start_eden():
    eden = {(0, 0, 0): [1, 0, 0],
            (0, 0, 1): [0, 1, 0],
            (0, 0, -1): [0, 1, 0],
            (1, 0, 0): [0, 1, 0],
            (-1, 0, 0): [0, 1, 0],
            (0, 1, 0): [0, 1, 0],
            (0, -1, 0): [0, 1, 0]}
    perimeter = [(0, 0, 1), (0, 0, -1), (1, 0, 0), (-1, 0, 0), (0, 1, 0), (0, -1, 0)]
    return eden, perimeter

def actualize_neighbors(tile_selected, eden, perimeter):
    n3 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2]]
    n4 = [tile_selected[0] + -1, tile_selected[1], tile_selected[2]]
    n5 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2]]
    n6 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2]]
    n1 = [tile_selected[0], tile_selected[1], tile_selected[2] + 1]
    n2 = [tile_selected[0], tile_selected[1], tile_selected[2] - 1]

    n1 = tuple(n1)
    n2 = tuple(n2)
    n3 = tuple(n3)
    n4 = tuple(n4)
    n5 = tuple(n5)
    n6 = tuple(n6)

    nearest_n_tiles = [n1, n2, n3, n4, n5, n6]

    nearest_n = [0, 0, 0, 0, 0, 0]

    if n1 in eden:
        eden[n1][1] = eden[n1][1] + 1
        if eden[n1][0] == 1:
            nearest_n[0] = 1
    else:
        eden[n1] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n1]

    if n2 in eden:
        eden[n2][1] = eden[n2][1] + 1
        if eden[n2][0] == 1:
            nearest_n[1] = 1
    else:
        eden[n2] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n2]

    if n3 in eden:
        eden[n3][1] = eden[n3][1] + 1
        if eden[n3][0] == 1:
            nearest_n[2] = 1
    else:
        eden[n3] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n3]

    if n4 in eden:
        eden[n4][1] = eden[n4][1] + 1
        if eden[n4][0] == 1:
            nearest_n[3] = 1
    else:
        eden[n4] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n4]

    if n5 in eden:
        eden[n5][1] = eden[n5][1] + 1
        if eden[n5][0] == 1:
            nearest_n[4] = 1
    else:
        eden[n5] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n5]

    if n6 in eden:
        eden[n6][1] = eden[n6][1] + 1
        if eden[n6][0] == 1:
            nearest_n[5] = 1
    else:
        eden[n6] = [0, 1, eden[tile_selected][2]]
        perimeter = perimeter + [n6]

    return eden, perimeter, nearest_n, nearest_n_tiles

def neighbours(eden, tile_selected):
    ###### For the Euler Characteristic

    #### Level one z = z -1, Level two z = 0, Level 3 Z = +1

    l11 = [tile_selected[0] - 1, tile_selected[1] + 1, tile_selected[2] - 1]
    l16 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2] - 1]
    l12 = [tile_selected[0] + 1, tile_selected[1] + 1, tile_selected[2] - 1]
    l15 = [tile_selected[0] - 1, tile_selected[1], tile_selected[2] - 1]
    l17 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2] - 1]
    l13 = [tile_selected[0] - 1, tile_selected[1] - 1, tile_selected[2] - 1]
    l18 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2] - 1]
    l14 = [tile_selected[0] + 1, tile_selected[1] - 1, tile_selected[2] - 1]

    l21 = [tile_selected[0] - 1, tile_selected[1] + 1, tile_selected[2]]
    l22 = [tile_selected[0] + 1, tile_selected[1] + 1, tile_selected[2]]
    l23 = [tile_selected[0] - 1, tile_selected[1] - 1, tile_selected[2]]
    l24 = [tile_selected[0] + 1, tile_selected[1] - 1, tile_selected[2]]

    l31 = [tile_selected[0] - 1, tile_selected[1] + 1, tile_selected[2] + 1]
    l36 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2] + 1]
    l32 = [tile_selected[0] + 1, tile_selected[1] + 1, tile_selected[2] + 1]
    l35 = [tile_selected[0] - 1, tile_selected[1], tile_selected[2] + 1]
    l37 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2] + 1]
    l33 = [tile_selected[0] - 1, tile_selected[1] - 1, tile_selected[2] + 1]
    l38 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2] + 1]
    l34 = [tile_selected[0] + 1, tile_selected[1] - 1, tile_selected[2] + 1]

    nn = [tuple(l11), tuple(l16), tuple(l12), tuple(l15), tuple(l17), tuple(l13), tuple(l18), tuple(l14),
          tuple(l21), tuple(l22), tuple(l23), tuple(l24), tuple(l31), tuple(l36), tuple(l32), tuple(l35),
          tuple(l37), tuple(l33), tuple(l38), tuple(l34)]
    n = [0] * 20
    for i in range(0, len(nn)):
        if nn[i] in eden:
            if eden[nn[i]][0] == 1:
                n[i] = 1
    return n

def actualize_vef(vertices, edges, faces, nearest_n, n):
    v = [1] * 8
    e = [1] * 12
    f = [1] * 6

    if n[0] == 1:
        v[0] = 0
    if n[1] == 1:
        v[0] = 0
        v[1] = 0
        e[0] = 0
    if n[2] == 1:
        v[1] = 0
    if n[3] == 1:
        v[0] = 0
        v[2] = 0
        e[3] = 0
    if n[4] == 1:
        v[1] = 0
        v[3] = 0
        e[1] = 0
    if n[5] == 1:
        v[2] = 0
    if n[6] == 1:
        v[2] = 0
        v[3] = 0
        e[2] = 0
    if n[7] == 1:
        v[3] = 0

    if n[8] == 1:
        v[0] = 0
        v[4] = 0
        e[11] = 0
    if n[9] == 1:
        v[1] = 0
        v[5] = 0
        e[8] = 0
    if n[10] == 1:
        v[2] = 0
        v[6] = 0
        e[10] = 0
    if n[11] == 1:
        v[3] = 0
        v[7] = 0
        e[9] = 0

    if n[12] == 1:
        v[4] = 0
    if n[13] == 1:
        v[4] = 0
        v[5] = 0
        e[4] = 0
    if n[14] == 1:
        v[5] = 0
    if n[15] == 1:
        v[4] = 0
        v[6] = 0
        e[7] = 0
    if n[16] == 1:
        v[5] = 0
        v[7] = 0
        e[5] = 0
    if n[17] == 1:
        v[6] = 0
    if n[18] == 1:
        v[6] = 0
        v[7] = 0
        e[6] = 0
    if n[19] == 1:
        v[7] = 0

    if nearest_n[0] == 1:
        v[4] = 0
        v[5] = 0
        v[7] = 0
        v[6] = 0
        e[4] = 0
        e[5] = 0
        e[6] = 0
        e[7] = 0
        f[0] = 0
    if nearest_n[1] == 1:
        v[0] = 0
        v[1] = 0
        v[2] = 0
        v[3] = 0
        e[0] = 0
        e[1] = 0
        e[2] = 0
        e[3] = 0
        f[1] = 0
    if nearest_n[2] == 1:
        v[1] = 0
        v[5] = 0
        v[3] = 0
        v[7] = 0
        e[1] = 0
        e[5] = 0
        e[8] = 0
        e[9] = 0
        f[2] = 0
    if nearest_n[3] == 1:
        v[0] = 0
        v[4] = 0
        v[2] = 0
        v[6] = 0
        e[3] = 0
        e[7] = 0
        e[11] = 0
        e[10] = 0
        f[3] = 0
    if nearest_n[4] == 1:
        v[0] = 0
        v[4] = 0
        v[1] = 0
        v[5] = 0
        e[0] = 0
        e[4] = 0
        e[11] = 0
        e[8] = 0
        f[4] = 0
    if nearest_n[5] == 1:
        v[2] = 0
        v[6] = 0
        v[3] = 0
        v[7] = 0
        e[2] = 0
        e[6] = 0
        e[9] = 0
        e[10] = 0
        f[5] = 0

    vertices = vertices + sum(v)
    edges = edges + sum(e)
    faces = faces + sum(f)

    return vertices, edges, faces, sum(v), sum(e), sum(f)

def euler_characteristic(k0, k1, k2, k3):
    return k0 - k1 + k2 - k3

def add_neighbours_bds(bds, j, iterations, num_possible_components, merged, finished, eden):
    tile_selected = bds[j][iterations]

    n3 = [tile_selected[0] + 1, tile_selected[1], tile_selected[2]]
    n4 = [tile_selected[0] + -1, tile_selected[1], tile_selected[2]]
    n5 = [tile_selected[0], tile_selected[1] + 1, tile_selected[2]]
    n6 = [tile_selected[0], tile_selected[1] - 1, tile_selected[2]]
    n1 = [tile_selected[0], tile_selected[1], tile_selected[2] + 1]
    n2 = [tile_selected[0], tile_selected[1], tile_selected[2] - 1]

    n1 = tuple(n1)
    n2 = tuple(n2)
    n3 = tuple(n3)
    n4 = tuple(n4)
    n5 = tuple(n5)
    n6 = tuple(n6)

    nearest_n_tiles = [n1, n2, n3, n4, n5, n6]

    nearest_n = [0, 0, 0, 0, 0, 0]

    if n1 in eden:
        if eden[n1][0] == 0:
            nearest_n[0] = 1
    else:
        nearest_n[0] = 1

    if n2 in eden:
        if eden[n2][0] == 0:
            nearest_n[1] = 1
    else:
        nearest_n[1] = 1

    if n3 in eden:
        if eden[n3][0] == 0:
            nearest_n[2] = 1
    else:
        nearest_n[2] = 1

    if n4 in eden:
        if eden[n4][0] == 0:
            nearest_n[3] = 1
    else:
        nearest_n[3] = 1

    if n5 in eden:
        if eden[n5][0] == 0:
            nearest_n[4] = 1
    else:
        nearest_n[4] = 1

    if n6 in eden:
        if eden[n6][0] == 0:
            nearest_n[5] = 1
    else:
        nearest_n[5] = 1

    for i in range(0, 6):
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

def increment_betti_2(eden, tile_selected, nearest_n, nearest_n_tiles, barcode, time, holes, total_holes, created_holes,
                      tags):

    total_holes_old = total_holes
    barcode_old = barcode
    holes_old = holes
    created_holes_old = created_holes
    tags_old = tags

    if eden[tile_selected][2] == 0:
        per = 1  # This is 1 if the tile added was in the out perimeter
    else:
        num_hole = eden[tile_selected][2]
        per = 0
    # In this case the tile added was in a hole

    betti_2 = 0

    if sum(nearest_n) == 6:
        betti_2 = - 1
        barcode[num_hole][1][1] = float(time + 2)  # Are we covering a hole that was never divided?
        holes[num_hole].remove(tile_selected)
        if holes[num_hole] == []:
            holes.pop(num_hole)

    if sum(nearest_n) == 5:
        betti_2 = 0
        if per == 0:
            holes[num_hole].remove(tile_selected)
    # print(nearest_n)
    if sum(nearest_n) != 6 and sum(nearest_n) != 5:
        num_possible_components = 0
        bds = []
        iterations = 0
        for i in range(0, 6):
            if nearest_n[i] == 0:
                num_possible_components = num_possible_components + 1
                bds = bds + [[nearest_n_tiles[i]]]

        finished = [0] * num_possible_components
        merged = finished.copy()

        while sum(finished) < num_possible_components - per:
            for j in range(0, num_possible_components):
                if finished[j] == 0:
                    bds, merged, finished = add_neighbours_bds(bds, j, iterations, num_possible_components, merged,
                                                               finished, eden)
                    if (iterations + 1) == len(bds[j]):
                        finished[j] = 1
            iterations = iterations + 1

        betti_2 = (num_possible_components - 1) - sum(merged)

        # print(betti_2, per)
        # At this point we have the bds components and the ones that were not merged will become the holes.
        # Here we actualize the holes and we actualize Hole No in eden.
        if betti_2 == 0:
            if per == 0:
                holes[num_hole].remove(tile_selected)

        else:
            if per == 1:
                for i in range(0, num_possible_components):
                    if finished[i] == 1 and merged[i] == 0:
                        total_holes = total_holes + 1
                        holes[total_holes] = bds[i].copy()

                        for x in bds[i]:
                            if x in eden:
                                eden[x][2] = total_holes

                        barcode[total_holes] = [0, [float(time + 2), float(0)], [total_holes]]
                        created_holes = created_holes + [[barcode[total_holes][2], bds[i].copy(), len(bds[i])]]

            else:
                if barcode[num_hole][0] == 0:
                    tags = tags + [num_hole]
                    barcode[num_hole][0] = 1

                holes.pop(num_hole)

                for i in range(0, num_possible_components):
                    if finished[i] == 1 and merged[i] == 0:
                        total_holes = total_holes + 1
                        holes[total_holes] = bds[i].copy()
                        for x in bds[i]:
                            if x in eden:
                                eden[x][2] = total_holes
                        barcode[total_holes] = [1, [float(time + 2), float(0)], barcode[num_hole][2] + [total_holes]]
                        created_holes = created_holes + [[barcode[total_holes][2], bds[i].copy(), len(bds[i])]]

    return betti_2, total_holes, eden, barcode, holes, created_holes, tags

def return_betti_1(betti_2, euler_ch):
    # X = V - E + F - Cubes
    # X = b_0 - b_1 + b_2
    # Observe that b_0 is 1 because polyominoes are connected
    # Cubes us exactly the time because we add one cube at each time
    return 1 + betti_2 - euler_ch


##########################
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
        leaves_parent = [x for x in b if len(x) == n - 1]

        #        print(leaves_parent)
        possible_leaves = [x for x in b if len(x) == n]
        #        print(possible_leaves)
        for j in leaves_parent:
            leaves = []
            for x in possible_leaves:
                #            print(x)
                root = list(x)
                del root[-1]
                root = tuple(root)
                if hamming2(j, root) == 0:
                    leaves = leaves + [x]
            #            Diff(possible_leaves, leaves)
            #            print(leaves)
            if len(leaves) > 0:
                times = []
                for x in leaves:
                    times = times + [b[x][1]]
                if 0 in times:
                    ind = times.index(0)
                    for i in range(0, len(leaves)):
                        if i != ind:
                            bars = bars + [b[leaves[i]]]
                else:
                    ind = times.index(max(times))
                    for i in range(0, len(leaves)):
                        if i != ind:
                            bars = bars + [b[leaves[i]]]
                    b[j][1] = max(times)
        n = n - 1
    bars = bars + [b[(tag,)]]
    return bars



# Time = 10000
# Eden, Perimeter, Betti_2_total_vector, Betti_2_vector_changes, Barcode, Holes, Betti_1_total, \
#     Betti_1_total_vector, Created_holes, Process, Perimeter_len, Skipped, I, Final_barcode = grow_eden(Time)
#
# convert_perseus(Process)
#
# eden_model = gd.CubicalComplex(perseus_file='300000_3D_1_final.txt')
#
# eden_model.persistence()
# A = eden_model.persistence_intervals_in_dimension(1)
# B = [elem for elem in A if elem[1] == float('inf')]
# A = eden_model.persistence_intervals_in_dimension(2)
# B = [elem for elem in A if elem[1] == float('inf')]
# final = np.array(Final_barcode)
# A_sorted = A.sort()
# final_sorted = final.sort()
# print(A_sorted == final_sorted)
# # gd.plot_persistence_diagram(eden_model.persistence(), legend=True)
# # gd.plot_persistence_barcode(eden_model.persistence())
# gd.plot_persistence_density(eden_model.persistence(), legend=True)
# a = 10
