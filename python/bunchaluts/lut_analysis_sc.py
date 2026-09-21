# author: Suhas Vittal

from sys import argv
import math as m

d = int(argv[1])
C = float(argv[2])

R = d

# 1. Compute amount of detection events in window
N = (d*d-1)//2 * R
print(f'Events = {N}')

# 2. Assuming physical error rate (`p`) is 1e-3, compute probability and number of 1-errors, 2-errors, etc.
#   cutoff at `C`
p = 1e-4
num_events = 0
Pr = 0.0

k = 0
while Pr < 1-C and k < 100:
    cnt = m.comb(N,k)
    Pr += cnt * (1-p)**(N-k) * p**k
    num_events = cnt
    k += 1

print(f'Threshold k < {k}, Pr = {Pr}, events = {num_events}')
