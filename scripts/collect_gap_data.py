import os

SHOTS = 100_000_000
MAX_ERRORS = 10_000
CORES = 80
P = 0.0015

for d in [7, 9, 11, 13]:
    cmd = f'mpirun -np {CORES} ./build/qudec pymatching {d} -r {3*d} -p {P} -m monte_carlo --mc-max-samples {SHOTS} --mc-stop-at-errors {MAX_ERRORS} -pp -g > out/gap/d{d}.out'
    print(cmd)
    os.system(cmd)

cmd = f'mpirun -np {CORES} ./build/qudec pymatching 23 -r {23*3} -p {P} -m monte_carlo --mc-max-samples 100000 -pp -g > out/gap/d23_distr_only.out'
print(cmd)
os.system(cmd)

