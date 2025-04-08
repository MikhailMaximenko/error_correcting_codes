#!/bin/bash

rm -f RM_8_4_0.txt RM_8_4_1.txt RM_8_4_2.txt RM_32_16_0.txt RM_32_16_1.txt RM_32_16_2.txt

./build/task-exec RM_8_4.gen RM_8_4_0.txt 0 100000 0.3
./build/task-exec RM_8_4.gen RM_8_4_1.txt 1 100000 0.3
./build/task-exec RM_8_4.gen RM_8_4_2.txt 2 100000 0.3
./build/task-exec RM_32_16.gen RM_32_16_0.txt 0 100000 0.3
./build/task-exec RM_32_16.gen RM_32_16_1.txt 1 100000 0.3
./build/task-exec RM_32_16.gen RM_32_16_2.txt 2 100000 0.3

python3 build_plots.py << EOF
	RM_8_4_0.txt RM_8_4_1.txt RM_8_4_2.txt RM_32_16_0.txt RM_32_16_1.txt RM_32_16_2.txt
EOF
