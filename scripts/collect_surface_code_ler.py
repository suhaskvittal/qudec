import os

DISTANCES = [5,7,9,11,13]
P = [3e-3,2e-3,1.5e-3,1e-3]
CORES = 64
DECODERS = [  # output folder prefix, decoder name, options
    ('uf', 'uf', ''),
    ('pymatching', 'pymatching', ''),
    ('correlated_pymatching', 'pymatching', '--correlated'),
    ('libra', 'seeded', '--seed-decoder pymatching --second-decoder libra --seed-gap cluster --seed-gap-threshold-db 40'),
    ('seeded_correlated_pymatching', 'seeded', '--seed-decoder uf --second-decoder correlated_pymatching')
]

SHOTS = 10_000_000
COMMON_OPTS = f'--mc-max-samples {SHOTS} -pp'

for (folder_prefix, decoder_name, opts) in DECODERS:
    output_folder = f'out/{folder_prefix}_ler'
    os.system(f'mkdir -p {output_folder}')
    for d in DISTANCES:
        for p in P:
            output_file = f'{output_folder}/d{d}_p{p}.out'
            cmd = f'mpirun -np {CORES} ./build/qudec {decoder_name} -d {d} -r {d} {COMMON_OPTS} {opts} > {output_file}'
            print(cmd)
            os.system(cmd)
