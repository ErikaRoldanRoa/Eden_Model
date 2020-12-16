#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import random
from tqdm import tqdm
import os
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import collections
from scipy.optimize import curve_fit
import matplotlib.patches as mpatches
import gudhi as gd


"""
Created on Thu Feb 13 19:04:05 2020
@author: Erika Roldan Roa
"""
"""
Created on Thu Nov 21 19:34:33 2019
@author: err
"""

"""GROWING"""
def grow_eden(t):
    vertices = 4
    edges = 4

    eden, perimeter = start_eden_2d()  # perimeter is an array consisting of all tiles that are on the perimeter
    process = [(0, 0)]  # an array consisting of all tiles that were added

    perimeter_len = []  # an array consisting of perimeter lengths on every time step
    holes = {}
    total_holes = 0
    barcode = {}
    created_holes = []
    tags = []

    betti_1_total = 0
    betti_1_vector_changes = [0]
    betti_1_total_vector = [0]
    pbar = tqdm(total=t, position=0, leave=True)
    pbar.update(1)
    for i in (range(1, t)):
        pbar.update(1)

        perimeter_len = perimeter_len + [len(perimeter)]
        x = random.randint(0, len(perimeter) - 1)

        tile_selected = perimeter[x]
        process = process + [tile_selected]

        perimeter.pop(x)
        eden[tile_selected][0] = 1

        eden, perimeter, nearest_n, nearest_neighbour_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        nearest_diag_tiles = neighbours_diag(eden, tile_selected)

        vertices, edges = actualize_vef(vertices, edges, nearest_n, nearest_diag_tiles)
        betti_1, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_1(eden, tile_selected,
                                                                                            nearest_n,
                                                                                            nearest_neighbour_tiles, barcode, i,
                                                                                            holes, total_holes,
                                                                                            created_holes, tags)
        # print('betti_1: ', betti_1)
        # draw_polyomino(eden, 0)
        betti_1_vector_changes = betti_1_vector_changes + [betti_1]
        betti_1_total += betti_1
        betti_1_total_vector = betti_1_total_vector + [betti_1_total]

    pbar.close()
    final_barcode = barcode_forest(barcode, tags)

    l = len(perimeter)
    perimeter_len = perimeter_len + [l]

    return eden, perimeter, betti_1_total_vector, betti_1_vector_changes, barcode, holes, betti_1_total, \
           created_holes, process, perimeter_len, tags, final_barcode

def grow_eden_debugging(t, ordered_tiles):
    vertices = 4
    edges = 4

    eden, perimeter = start_eden_2d()

    holes = {}
    total_holes = 0
    barcode = {}
    created_holes = []
    tags = []

    betti_1_total = 0
    betti_1_vector = []
    betti_1_euler_total = 0
    betti_1_total_vector = [0]
    len_perimeter = [4]
    pbar = tqdm(total=t, position=0, leave=True)
    pbar.update(1)
    for i in (range(1, t)):
        pbar.update(1)

        tile_selected = ordered_tiles[i]
        perimeter.remove(tile_selected)
        eden[tile_selected][0] = 1

        eden, perimeter, nearest_n, nearest_n_tiles = actualize_neighbors(tile_selected, eden, perimeter)
        n = neighbours_diag(eden, tile_selected)

        vertices, edges = actualize_vef(vertices, edges, nearest_n, n)
        betti_1, total_holes, eden, barcode, holes, created_holes, tags = increment_betti_1(eden, tile_selected,
                                                                                            nearest_n,
                                                                                            nearest_n_tiles, barcode, i,
                                                                                            holes, total_holes,
                                                                                            created_holes, tags)
        betti_1_vector = betti_1_vector + [betti_1]
        betti_1_total = betti_1_total + betti_1
        betti_1_total_vector += [betti_1_total]

        betti_1_euler_total = increment_betti_1_euler(vertices, edges, i)
        if betti_1_total != betti_1_euler_total:
            raise ValueError('betti_1_total does not equal betti_1 with euler')
        len_perimeter += [len(perimeter)]

    pbar.close()
    final_barcode = barcode_forest(barcode, tags)

    return eden, perimeter, betti_1_vector, betti_1_total_vector, barcode, holes, betti_1_total, betti_1_euler_total, created_holes, tags, final_barcode, len_perimeter

"""GUDHI"""
def convert_gudhi(process, folder_name):
    min_x = min(process, key=lambda x: x[0])
    max_x = max(process, key=lambda x: x[0])
    min_y = min(process, key=lambda x: x[1])
    max_y = max(process, key=lambda x: x[1])

    dimension = 2

    long = max_x[0] + ((-1) * (min_x[0])) + 1
    wide = max_y[1] + ((-1) * (min_y[1])) + 1

    time = len(process)
    filename = folder_name+'/gudhi_'+str(time)+'.txt'

    total = long*wide
    pbar = tqdm(total=total)

    with open(filename, 'w+') as f:
        f.writelines('%d\n' % dimension)
        f.writelines('%d\n' % long)
        f.writelines('%d\n' % wide)

        for i in range(min_y[1], max_y[1] + 1):
            for j in range(min_x[0], max_x[0] + 1):
                pbar.update(1)
                if (j, i) in process:
                    f.writelines('%d\n' % process.index((j, i)))
                else:
                    f.writelines('inf\n')
    pbar.close()
    return filename

def gudhi_analysis(filename, final_barcode, folder_name, length):
    print("Creating Cubical Complex...")
    eden_model = gd.CubicalComplex(perseus_file=filename)
    eden_model.persistence()
    A = eden_model.persistence_intervals_in_dimension(1)
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
    plt.close()

"""DRAWING"""
def draw_square(x, y, col='gray', ls=0.35):
    """With center at x, y draw a square of area 1"""
    """it's area is actually 4ls^2, isn't it?"""
    """ls is a half of square side"""
    plt.fill([x - ls, x + ls, x + ls, x - ls], [y - ls, y - ls, y + ls, y + ls], color=col)

def draw_polyomino(eden, time, folder_name):
    """This function draws a square for each (x, y) if the respective IO entry is 1 and does nothing if it is 0"""
    """As I got entries of eden are (x, y): [a, b, c] and if a == 1 than we draw the square with the center (x, y)
    b probably corresponds to the number how many times we could've add this square to the eden"""
    plt.style.use('ggplot')
    plt.axis('off')
    plt.rc('grid', linestyle="-", color='black')
    plt.grid(True)
    # plt.rc('grid', linestyle="-", color='black')  # why do we need this line twice? and does it change anything?
    plt.gca().set_aspect('equal', adjustable='box')
    pbar = tqdm(total=time, position=0, leave=True)
    for x in eden:
        if eden[x][0] == 1:
            pbar.update(1)
            draw_square(x[0], x[1], 'gray')
    draw_square(0, 0, 'green')
    plt.savefig(folder_name+'/eden.svg', format='svg', dpi=1200)
    # plt.show()
    plt.close()

def draw_polyomino_holes(eden, holes, time):
    """This function draws a square for each (x, y) if the respective IO entry is 1 and does nothing if it is 0"""
    plt.style.use('ggplot')
    plt.axis('off')
    plt.grid(True)
    plt.gca().set_aspect('equal', adjustable='box')
    for x in eden:
        if eden[x][0] == 1:
            draw_square(x[0], x[1], 'gray')
    for x in holes:
        s = len(holes[x])
        for i in range(0, s):
            draw_square(holes[x][i][0], holes[x][i][1], 'red')
    draw_square(0, 0, 'green')
    plt.savefig('pictures/eden_' + str(time) + '_holes.svg', format='svg', dpi=1200)
    plt.show()

def draw_barcode(barcode, time):
    fig = plt.figure()
    plt.style.use('ggplot')
    # plt.axis('off')
    # plt.grid(True)
    plt.rc('grid', linestyle="-", color='gray')
    plt.yticks([])
    # plt.gca().set_aspect('equal', adjustable='box')
    i = 0
    for x in barcode:
        if barcode[x][1][1] == 0:
            plt.plot([barcode[x][1][0], time], [i, i], 'k-', lw=2)
        else:
            plt.plot([barcode[x][1][0], barcode[x][1][1]], [i, i], 'k-', lw=2)
        i = i + 40
    fig.savefig('2d/'+str(int(time/1000))+'k/barcode_'+str(time)+'.svg', format='svg', dpi=1200)
    plt.show()

def draw_frequencies_1(dict, changes, folder_name):
    fig = plt.figure()
    ax = fig.add_subplot(1, 1, 1)
    l = len(dict[0])

    ch_4 = [i for i, j in enumerate(changes) if j == 4]
    y_4 = []
    for x in ch_4:
        y_4 += [dict[4][x+1]]

    sh = []
    for j in np.arange(-1, 3):
        sh.append(next((i for i, x in enumerate(dict[j]) if x), 0))
    shift = max(sh)

    ax.plot(range(shift, l), dict[-1][shift:], color='tab:red', label='-1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[0][shift:], color='tab:orange', label='0',  linewidth=0.75)
    ax.plot(range(shift, l), dict[1][shift:], color='tab:green', label='+1',  linewidth=0.75)
    ax.plot(range(shift, l), dict[2][shift:], color='tab:blue', label='+2',  linewidth=0.75)

    plt.yscale('log')
    # ax.set_title('betti_1 frequencies')
    ax.set_ylabel(r'Frequency of Change in $\beta_1$')
    ax.set_xlabel('t')
    ax.ticklabel_format(style='sci', axis='x', scilimits=(0, 0))
    # ax.ticklabel_format(useOffset=False)
    ax.legend(loc=1, prop={'size': 6})
    fig.savefig(folder_name+'/fr_b_1_.png', format='png', dpi=1200)
    # plt.show()
    plt.close()

"""SUPPLEMENTARY FUNCTIONS"""
def return_frequencies_1(vect, time):
    changes = [vect[i+1]-vect[i] for i in range(len(vect)-1)]
    # values = list(set(changes))
    # values.sort()
    values = [-1, 0, 1, 2]
    freq = {i : [0] for i in values}

    for i in tqdm(range(1, time+1), position=0, leave=True):
        counter = collections.Counter(changes[:i])
        for k in values:
            freq[k].append(counter[k]/i)
    return freq, changes

def read_eden_txt(filename):
    eden = []
    for t in open(filename).read().split('), ('):
        # print(t)
        split = t.strip('()[]').split(',')
        split[-2] = split[-2].strip(')')
        # split[-1] = split[-1].strip(")]\n")
        tup = tuple([int(x) for x in split[:-1]])
        eden.append((tup, float(split[-1])))
    return eden

def read_eden_perseus(filename):
    eden = []
    file = open(filename, 'r')
    lines = file.readlines()
    # dim = int(lines[0])
    for line in lines[1:]:
        split = line.strip('\n').split(' ')
        tup = tuple([int(x) for x in split[:-1]])
        eden.append((tup, float(split[-1])))
    return eden

def hamming2(s1, s2):
    """Calculate the Hamming distance between two bit strings"""
    assert len(s1) == len(s2)
    return sum(c1 != c2 for c1, c2 in zip(s1, s2))

def diff(li1, li2):
    """Returns characters that present only in one of two strings"""
    """the function isn't used (the line with it is commented) do we need it?"""
    li_dif = [i for i in li1 + li2 if i not in li1 or i not in li2]
    return li_dif

def start_eden_2d():
    """eden is a dictionary consisted of such items: (x, y): [a, b, c]
    where (x, y) are the coordinates of the square center
    a = 1 if the square is already in the complex (0 if only in the perimeter)
    b = number of neighbours already in the complex (from 0 to 4)
    c = 1 if tile is a part of a hole"""
    """perimeter is a layer of the squares lying on the perimeter (but not yet it the complex)"""
    eden = {(0, 0): [1, 0, 0],
            (1, 0): [0, 1, 0],
            (-1, 0): [0, 1, 0],
            (0, 1): [0, 1, 0],
            (0, -1): [0, 1, 0]}
    perimeter = [(1, 0), (-1, 0), (0, 1), (0, -1)]
    return eden, perimeter

def actualize_neighbors(tile_selected, eden, perimeter):

    n2 = (tile_selected[0] + 1, tile_selected[1])
    n4 = (tile_selected[0] + -1, tile_selected[1])
    n1 = (tile_selected[0], tile_selected[1] + 1)
    n3 = (tile_selected[0], tile_selected[1] - 1)

    nearest_n_tiles = [n1, n2, n3, n4]  # coordinates of the neighbours
    nearest_n = [0, 0, 0, 0]  # 1 if neighbour in the complex, else: 0
    for i, n in enumerate(nearest_n_tiles):
        if n in eden:
            eden[n][1] = eden[n][1] + 1
            if eden[n][0] == 1:
                nearest_n[i] = 1
        else:
            eden[n] = [0, 1, eden[tile_selected][2]]
            perimeter = perimeter + [n]

    return eden, perimeter, nearest_n, nearest_n_tiles

def neighbours_diag(eden, tile_selected):
    """For the Euler Characteristic
    Level one z = z - 1, Level two z = 0, Level three Z = +1"""

    l24 = (tile_selected[0] - 1, tile_selected[1] + 1)
    l21 = (tile_selected[0] + 1, tile_selected[1] + 1)
    l23 = (tile_selected[0] - 1, tile_selected[1] - 1)
    l22 = (tile_selected[0] + 1, tile_selected[1] - 1)

    N = [l21, l22, l23, l24]
    nearest_diag = [0] * 4
    for i in range(0, len(N)):
        if N[i] in eden:
            if eden[N[i]][0] == 1:
                nearest_diag[i] = 1
    return nearest_diag

def actualize_vef(vertices, edges, nearest_n, nearest_diag):
    """for Euler characteristics. function updates number of vertices and edges"""
    # CHANGE! the function was wrong. it was calculating the number of vertices in a wrong way
    v = [1] * 4
    e = [1] * 4
    for i in range(4):
        if nearest_diag[i] == 1:
            v[i] = 0
        if nearest_n[i] == 1:
            e[i] = 0
            v[i-1] = 0
            v[i] = 0
    vertices = vertices + sum(v)
    edges = edges + sum(e)

    return vertices, edges

def add_neighbours_bds(bds, j, iterations, num_possible_components, merged, finished, eden):  # spreading gas
    """j is number of neighbour whose gas we are spreading"""
    tile_selected = bds[j][iterations]

    n1 = (tile_selected[0] + 1, tile_selected[1])
    n2 = (tile_selected[0] + -1, tile_selected[1])
    n3 = (tile_selected[0], tile_selected[1] + 1)
    n4 = (tile_selected[0], tile_selected[1] - 1)

    nearest_n_tiles = [n1, n2, n3, n4]
    nearest_n = [0, 0, 0, 0]  # equals 1 if corresponding neighbour is 'free', i.e. not in complex

    for i, n in enumerate(nearest_n_tiles):
        if n in eden:
            if eden[n][0] == 0:
                nearest_n[i] = 1
        else:
            nearest_n[i] = 1

    for i in range(0, 4):
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

def increment_betti_1(eden, tile_selected, nearest_n, nearest_n_tiles, barcode, time, holes, total_holes, created_holes,
                      tags):
    if eden[tile_selected][2] == 0:
        per = 1  # This is 1 if the tile added was in the out perimeter
    else:
        num_hole = eden[tile_selected][2]
        per = 0
    # In this case the tile added was in a hole
    betti_1 = 0

    if sum(nearest_n) == 4:
        betti_1 = - 1
        barcode[num_hole][1][1] = time + 1  # Are we covering a hole that was never divided?
        holes[num_hole].remove(tile_selected)
        if not holes[num_hole]:
            holes.pop(num_hole)

    if sum(nearest_n) == 3:
        betti_1 = 0
        if per == 0:
            holes[num_hole].remove(tile_selected)
    # print(nearest_n)
    if sum(nearest_n) < 3:
        num_possible_components = 0
        bds = []
        iterations = 0
        for i in range(0, 4):
            if nearest_n[i] == 0:
                num_possible_components += 1
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

        betti_1 = (num_possible_components - 1) - sum(merged)  # imp: it's not a betti_1 but a change in betti_1
        # print(betti_1, per)
        # At this point we have the bds components and the ones that were not merged will become the holes.
        # Here we actualize the holes and we actualize Hole No in eden.
        if betti_1 == 0:
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
                        barcode[total_holes] = [0, [time + 1, 0], [total_holes]]
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
                        if bds[i] == [(-36, 85), (-36, 84), (-35, 84), (-34, 84)]:
                            a = 10
                        for x in bds[i]:
                            if x in eden:
                                eden[x][2] = total_holes
                        barcode[total_holes] = [1, [time + 1, 0], barcode[num_hole][2] + [total_holes]]
                        created_holes = created_holes + [[barcode[total_holes][2], bds[i].copy(), len(bds[i])]]

    return betti_1, total_holes, eden, barcode, holes, created_holes, tags

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
        # print(b)
        bars_hole = bars_hole + bars_from_tree(b, x)
    return bars_pure + bars_hole

def bars_from_tree(b, tag):
    n = max(len(x) for x in b)
    bars = []
    while n > 0:
        leaves_parent = [x for x in b if len(x) == n - 1]

        # print(leaves_parent)
        possible_leaves = [x for x in b if len(x) == n]
        # print(possible_leaves)
        for j in leaves_parent:
            leaves = []
            for x in possible_leaves:
                # print(x)
                root = list(x)
                del root[-1]
                root = tuple(root)
                if hamming2(j, root) == 0:
                    leaves = leaves + [x]
            # diff(possible_leaves, leaves)
            # print(leaves)
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

def increment_betti_1_euler(vertices, edges, time):
    return 1 - vertices + edges - time - 1

def num_holes(created_holes, holes):
    tromino = dict.fromkeys(['l', 'i'], 0)
    tromino_f = dict.fromkeys(['l', 'i'], 0)
    tetromino = dict.fromkeys(['I', 'L', 'T', 'O', 'Z'], 0)
    tetromino_f = dict.fromkeys(['I', 'L', 'T', 'O', 'Z'], 0)

    total_3 = [i for i in created_holes if i[-1] == 3]
    final_3 = [holes[i] for i in holes if len(holes[i]) == 3]
    total_4 = [i for i in created_holes if i[-1] == 4]
    final_4 = [holes[i] for i in holes if len(holes[i]) == 4]

    # tromino case
    for x in total_3:
        j = 0
        for i in range(2):
            if x[-2][0][i] == x[-2][1][i] == x[-2][2][i]:
                j += 1
        long = j
        tromino['i'] += long
        tromino['l'] += (1 - long)
    for x in final_3:
        j = 0
        for i in range(2):
            if x[0][i] == x[1][i] == x[2][i]:
                j += 1
        long = j
        tromino_f['i'] += long
        tromino_f['l'] += (1 - long)

    # tetromino case
    for x in total_4:
        dist = distances(x[-2])
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
        tetromino[typ] += 1
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
        tetromino_f[typ] += 1

    return tromino, tromino_f, tetromino, tetromino_f

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

"""PLOTTING"""

def draw_diagram_holes(created_holes, holes, folder_name):
    fr_cr = [created_holes[i][-1] for i in range(len(created_holes))]
    fr_cr.sort()
    fr_final = [len(holes[i]) for i in holes]
    fr_final.sort()

    counter_cr = collections.Counter(fr_cr)
    for j in range(1, list(counter_cr.keys())[-1]):
        if j not in counter_cr:
            counter_cr[j] = 0

    counter_final = collections.Counter(fr_final)
    for i in counter_cr.keys():
        if i not in counter_final:
            counter_final[i] = 0
    width = 0.35

    labels = range(len(counter_cr.keys())+1)
    x = np.arange(len(labels))

    fig, ax = plt.subplots()
    plt.yscale('log')
    ax.bar(np.array(list(counter_cr.keys())) - width/2, counter_cr.values(), width, color=[(0.44, 0.57, 0.79)], label='Total')
    ax.bar(np.array(list(counter_final.keys())) + width/2, counter_final.values(), width, color=[(225/256, 151/256, 76/256)], label='Final')

    ax.set_ylabel('Frequency of Number of Holes')
    ax.set_xlabel('Volume')
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    if len(labels) >= 60:
        plt.setp(ax.get_xticklabels(), fontsize=3)
    elif len(labels) >= 50:
        plt.setp(ax.get_xticklabels(), fontsize=4)
    elif len(labels) >= 45:
        plt.setp(ax.get_xticklabels(), fontsize=5)
    elif len(labels) >= 40:
        plt.setp(ax.get_xticklabels(), fontsize=6)
    elif len(labels) >= 30:
        plt.setp(ax.get_xticklabels(), fontsize=7)
    elif len(labels) >= 20:
        plt.setp(ax.get_xticklabels(), fontsize=6)

    ax.legend()
    fig.tight_layout()
    fig.savefig(folder_name+'/holes.png', format='png', dpi=1200)
    plt.close()

def plot_b_per(b1, p2, time, folder_name, times=None):
    n = int(time/10)
    if times is None:
        times = np.arange(1, time+1)

    def func(x, a, b):
        return a * x ** b
    ydata_f = b1
    xdata_f = times
    ydata = ydata_f
    xdata = xdata_f
    plt.xscale('log')
    plt.yscale('log')
    plt.plot(xdata_f[n:], ydata_f[n:], 'm-', label=r'$\beta_1(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata_f[n:], func(xdata_f[n:], *popt), 'm--', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt), linewidth=0.75)

    ydata = p2
    xdata = times
    plt.plot(xdata[n:], ydata[n:], color='orange', linestyle='solid', label=r'$P(t)$ data',  linewidth=0.75)
    popt, pcov = curve_fit(func, xdata, ydata)
    plt.plot(xdata[n:], func(xdata[n:], *popt), color='orange', linestyle='dashed', label=r'fit: $y=%5.2f x^{%5.3f}$' % tuple(popt),  linewidth=0.75)

    plt.xlabel('t')
    plt.ylabel('data')
    plt.legend(prop={'size': 6})
    plt.tight_layout()
    plt.savefig(folder_name+'/per-b-time_'+str(time)+'.png', dpi=1200)
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
    if er == [0, 0, 0, 0]:
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(loc='upper right')
    else:
        if er == [1, 1, 1, 1]:
            handles = []
        else:
            handles, labels2 = ax.get_legend_handles_labels()
        if er[0] == 1:
            patch = mpatches.Patch(color='navy', label='Trominoes Total', linewidth=0.35)
            handles.append(patch)
        if er[1] == 1:
            patch = mpatches.Patch(color='royalblue', label='Trominoes Final', linewidth=0.35)
            handles.append(patch)
        if er[2] == 1:
            patch = mpatches.Patch(color='chocolate', label='Tetrominoes Total', linewidth=0.35)
            handles.append(patch)
        if er[3] == 1:
            patch = mpatches.Patch(color='orange', label='Tetrominoes Final', linewidth=0.35)
            handles.append(patch)
        ax.set_xticks(x)
        ax.set_xticklabels(labels)
        ax.legend(handles=handles, loc='upper right')

    ax.set_ylabel('Frequency of Number of Holes')
    ax.set_xlabel('Type of a Hole')
    fig.tight_layout()
    fig.savefig(folder_name+'/tro-tetro-minoes.png', format='png', dpi=1200)
    plt.close()
