import matplotlib.pyplot as plt
import numpy as np
import math

files_with_stats = input().split()

plt.figure()
plt.xlabel('relation signal/noise (dB)')
plt.ylabel('logarithm error probability')
markers = ['o', 'v', '^', '<', '>', 's', '1', '2', '3', '4', '8', 'p']

cnt = 0
for file_name in files_with_stats:
    file = open(file=file_name, mode='r')

    stats = list(map(lambda x: x.rstrip(), file.readlines()))
    color = list(map(lambda x: x / 256, list(np.random.choice(range(256), 3))))
    errors = []
    signal_noises = []
    for row in stats:
        w, iters, hits, signal_noise = map(float, row.split())
        error = math.log10((iters - hits) / iters + 1e-9)
        errors.append(error)
        signal_noises.append(signal_noise)
    plt.plot(signal_noises, errors, marker=markers[cnt], color=color, label=file_name)
    cnt += 1


plt.legend(loc="upper right")
plt.savefig('results.png')
plt.show()
