#!/bin/bash


# ./build/task-exec RM_8_4.gen RM_8_4_0.txt 0 10 0.5
# ./build/task-exec RM_8_4.gen RM_8_4_1.txt 1 200 0.5
# ./build/task-exec RM_8_4.gen RM_8_4_2.txt 2 200 0.5
# ./build/task-exec ResIV.txt 5 1 0 0
# ./build/task-exec ResGV.txt 5 1 1 0
./build/task-exec ResGU.txt 5 1 1 1
# ./build/task-exec RM_32_16.gen RM_32_16_1.txt 1 200 0.5
# ./build/task-exec RM_32_16.gen RM_32_16_2.txt 2 200 0.5

# gnuplot -c make_plots.gp RM_8_4_0.txt RM_8_4_1.txt RM_8_4_2.txt RM_32_16_0.txt RM_32_16_1.txt RM_32_16_2.txt
