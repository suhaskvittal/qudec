import os 

CORES = 64
P = 1e-3

#for d in [5,7,9,11,29]:
for dec in ['pymatching', 'uf']:
    os.system(f'mkdir -p out/{dec}_gap')
    for d in [29]:
        shots = 100_000_000 if d <= 11 else 1_000_000
        if dec == 'uf' and d > 11:
            continue
        cmd = f'mpirun -np {CORES} ./build/qudec {dec} -d {d} -r {3*d} -p {P} --mc-max-samples {shots} --cluster-gap > out/{dec}_gap/d{d}.out'
        print(cmd)
        os.system(cmd)
