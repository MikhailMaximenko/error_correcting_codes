#!/bin/bash

./build/task-exec RM_64_22.gen trellis_based_rmld_results/ResIV_RM_64_22.txt 50 1 2 0 0
./build/task-exec RM_64_42.gen trellis_based_rmld_results/ResIV_RM_64_42.txt 50 1 2 0 0
./build/task-exec RM_64_57.gen trellis_based_rmld_results/ResIV_RM_64_57.txt 50 1 2 0 0
./build/task-exec RM_64_22.gen trellis_based_rmld_results/ResGV_RM_64_22.txt 50 1 2 1 0
./build/task-exec RM_64_42.gen trellis_based_rmld_results/ResGV_RM_64_42.txt 50 1 2 1 0
./build/task-exec RM_64_57.gen trellis_based_rmld_results/ResGV_RM_64_57.txt 50 1 2 1 0
./build/task-exec RM_64_22.gen trellis_based_rmld_results/ResGU_RM_64_22.txt 100 1 2 1 1
./build/task-exec RM_64_42.gen trellis_based_rmld_results/ResGU_RM_64_42.txt 100 1 2 1 1
./build/task-exec RM_64_57.gen trellis_based_rmld_results/ResGU_RM_64_57.txt 100 1 2 1 1
