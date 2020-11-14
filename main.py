import os
import numpy as np


def read_value(arr):
    while True:
        try:
            x = int(input())
            if x not in arr:
                raise ValueError
            break
        except ValueError:
            print("Oops!  That was no valid number.  Try again...")
    return x

print('Welcome to EDEN Model!')

print('Please, enter the desirable dimension of your model (from 2 to 5): ')
# dim = read_value([2, 3, 4, 5])
dim = 4

print('Do you have a file with a model? \n 1 -- you have a file \n 0 -- you want to generate a new model \n')
# file = read_value([0, 1])
file = 1

if dim <= 3:
    print('Do you want a picture of your model? (with a large a model it can take time) \n 1 -- yes \n 0 -- no \n')
    # pic = read_value([0, 1])
    pic = 0

"""NO FILE CASE"""
if file == 0:
    print('How many tiles would you like in your model?')
    # while True:
    #     try:
    #         Time = int(input())
    #         break
    #     except ValueError:
    #         print("Oops!  That was no valid number.  Try again...")
    Time = 10000
    if dim == 2:
        from e_2d import grow_eden, plot_b_per, draw_diagram_holes, num_holes, draw_tri_tetra, draw_barcode, \
            draw_polyomino, return_frequencies_1, draw_frequencies_1
        print("\nBuilding a model...")
        Eden, Perimeter, Betti_1_total_vector, Betti_1_vector_changes, Barcode, Holes, Betti_1_total, \
            Created_holes, Process, Perimeter_len, Tags, Final_barcode = grow_eden(Time)
        if not os.path.exists('2d/'+str(int(Time/1000))+'k/'):
            os.makedirs('2d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_1...")
        freq, changes = return_frequencies_1(Betti_1_total_vector, Time)
        draw_frequencies_1(freq, Time, changes)
        Tromino, Tromino_f, Tetromino, Tetromino_f = num_holes(Created_holes, Holes)
        draw_tri_tetra(Tromino, Tromino_f, Tetromino, Tetromino_f, Time)
        plot_b_per(Betti_1_total_vector, Perimeter_len, Time)

        # draw_barcode(Barcode, Time)
        if pic == 1:
            print("\nDrawing the complex...")
            draw_polyomino(Eden, Time)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
    if dim == 3:
        from e_3d import grow_eden, return_frequencies_1, draw_frequencies_1, num_holes, draw_tri_tetra, plot_b_per,\
            return_frequencies_2, draw_frequencies_2, grow_eden_debugging
        from e_2d import draw_diagram_holes
        print("\nBuilding a model...")
        Eden, Perimeter, Betti_2_total_vector, Betti_2_vector_changes, Barcode, Holes, Betti_1_total, \
            Betti_1_total_vector, Created_holes, Process, Perimeter_len, Skipped, I, Final_barcode = grow_eden(Time)

        # f = open("3d/sample_time_list.txt", "w+")
        # Process_file = [(tuple(x), i) for i, x in enumerate(Process)]
        # f.write(str(Process_file))
        # f.close()

        if not os.path.exists('3d/'+str(int(Time/1000))+'k/'):
            os.makedirs('3d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_1...")
        freq, changes = return_frequencies_1(Betti_1_total_vector, Time)
        draw_frequencies_1(freq, Time, changes)
        print("\nCalculating frequencies of betti_2...")
        freq, changes = return_frequencies_2(Betti_2_total_vector, Time)
        draw_frequencies_2(freq, Time, changes)
        Tricube, Tricube_f, Tetracube, Tetracube_f = num_holes(Created_holes, Holes)
        draw_tri_tetra(Tricube, Tricube_f, Tetracube, Tetracube_f, Time)
        plot_b_per(Betti_1_total_vector, Betti_2_total_vector, Perimeter_len, Time, 0)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
        if pic == 1:
            a = 1
            f = open("3d/"+str(int(Time/1000))+"k/MAYA.txt", "w+")
            f.write("import maya.cmds as cmds \nimport math as m \n"
                    "import os,sys \nEden = " + str(Process)+"\nt = len(Eden)"
                    "\nfor i in range(0,t):\n\taux = cmds.polyCube()"
                    "\n\tcmds.move(Eden[i][0],Eden[i][1],Eden[i][2],aux)")
            f.close()
            print("We created txt file \"MAYA\" for you. Just copy paste its content to MAYA!")
    if dim == 4:
        from e_4d import grow_eden, draw_frequencies_3, return_frequencies_3, plot_b_per
        from e_2d import draw_diagram_holes
        Eden, Perimeter, betti_3_vector, barcode, Holes, betti_3_total, Created_holes, Process, Perimeter_len,\
            Betti_3_total_vector = grow_eden(Time)
        # f = open("4d/sample_time_list.txt", "w+")
        # Process_file = [(tuple(x), i) for i, x in enumerate(Process)]
        # f.write(str(Process_file))
        # f.close()
        if not os.path.exists('4d/'+str(int(Time/1000))+'k/'):
            os.makedirs('4d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_3...")
        freq, changes = return_frequencies_3(Betti_3_total_vector, Time)
        draw_frequencies_3(freq, Time, changes)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
        plot_b_per(Betti_3_total_vector, Perimeter_len, Time, 0)


"""FILE CASE"""
if file == 1:
    if dim == 2:
        from e_2d import grow_eden_debugging, plot_b_per, draw_diagram_holes, num_holes, draw_tri_tetra, draw_barcode, draw_polyomino, \
            read_eden_txt, return_frequencies_1, draw_frequencies_1

        Eden_f = read_eden_txt("2d/sample_time_list.txt")
        Eden = [x[0] for x in Eden_f]
        Times = [x[1] for x in Eden_f]
        Time = len(Eden)
        print("\nComputing persistent homology...")
        Eden, Perimeter, Betti_1_vector, Betti_1_total_vector, Barcode, Holes, Betti_1_total, Betti_1_euler_total, \
            Created_holes, Tags, Final_barcode, Perimeter_len = grow_eden_debugging(len(Eden), Eden)
        if not os.path.exists('2d/'+str(int(Time/1000))+'k/'):
            os.makedirs('2d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_1...")
        freq, changes = return_frequencies_1(Betti_1_total_vector, Time)
        draw_frequencies_1(freq, Time, changes)
        Tromino, Tromino_f, Tetromino, Tetromino_f = num_holes(Created_holes, Holes)
        draw_tri_tetra(Tromino, Tromino_f, Tetromino, Tetromino_f, Time)
        plot_b_per(Betti_1_total_vector, Perimeter_len, Time, Times)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
        # draw_barcode(Barcode, Time)
        if pic == 1:
            print("Drawing the complex...")
            draw_polyomino(Eden, Time)
    if dim == 3:
        from e_3d import return_frequencies_1, draw_frequencies_1, num_holes, draw_tri_tetra, plot_b_per,\
            return_frequencies_2, draw_frequencies_2, grow_eden_debugging, read_eden_txt
        from e_2d import draw_diagram_holes
        Eden_f = read_eden_txt("3d/sample_time_list.txt")
        Eden = [x[0] for x in Eden_f]
        Times = [x[1] for x in Eden_f]
        Time = len(Eden)

        print("\nComputing persistent homology...")
        Eden, Perimeter, Betti_2_total_vector, Betti_1_total_vector, Barcode, Holes, \
            Betti_2_total, Betti_1_total, Created_holes, Perimeter_len = grow_eden_debugging(len(Eden), Eden)
        if not os.path.exists('3d/'+str(int(Time/1000))+'k/'):
            os.makedirs('3d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_1...")
        freq, changes = return_frequencies_1(Betti_1_total_vector, Time)
        draw_frequencies_1(freq, Time, changes)
        print("\nCalculating frequencies of betti_2...")
        freq, changes = return_frequencies_2(Betti_2_total_vector, Time)
        draw_frequencies_2(freq, Time, changes)
        Tricube, Tricube_f, Tetracube, Tetracube_f = num_holes(Created_holes, Holes)
        draw_tri_tetra(Tricube, Tricube_f, Tetracube, Tetracube_f, Time)
        plot_b_per(Betti_1_total_vector, Betti_2_total_vector, Perimeter_len, Time, 0)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
        if pic == 1:
            a = 1
            f = open("3d/"+str(int(Time/1000))+"k/MAYA.txt", "w+")
            f.write("import maya.cmds as cmds \nimport math as m \n"
                    "import os,sys \nEden = " + str(Process)+"\nt = len(Eden)"
                    "\nfor i in range(0,t):\n\taux = cmds.polyCube()"
                    "\n\tcmds.move(Eden[i][0],Eden[i][1],Eden[i][2],aux)")
            f.close()
            print("We created txt file \"MAYA\" for you. Just copy paste its content to MAYA!")
    if dim == 4:
        from e_4d import grow_eden_debugging, draw_frequencies_3, return_frequencies_3, plot_b_per, read_eden_txt
        from e_2d import draw_diagram_holes
        Eden_f = read_eden_txt("4d/sample_time_list.txt")
        Eden = [x[0] for x in Eden_f]
        Times = [x[1] for x in Eden_f]
        Time = len(Eden)
        print("\nComputing persistent homology...")
        Eden, Perimeter, betti_3_vector, barcode, Holes, betti_3_total, Created_holes, Perimeter_len,\
            Betti_3_total_vector = grow_eden_debugging(len(Eden), Eden)
        if not os.path.exists('4d/'+str(int(Time/1000))+'k/'):
            os.makedirs('4d/'+str(int(Time/1000))+'k/')
        print("\nCalculating frequencies of betti_3...")
        freq, changes = return_frequencies_3(Betti_3_total_vector, Time)
        draw_frequencies_3(freq, Time, changes)
        draw_diagram_holes(Created_holes, Holes, Time, dim)
        plot_b_per(Betti_3_total_vector, Perimeter_len, Time, 0)


print("WE ARE DONE!")
