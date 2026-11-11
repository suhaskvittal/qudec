import os

CORES = 16
DECODERS = [  # output folder prefix, decoder name, options
#   ('pymatching', 'pymatching', ''),
#   ('pymatching_with_gap', 'pymatching', '--cluster-gap'),
#   ('correlated_pymatching', 'pymatching', '--correlated'),
    ('seeded_correlated_pymatching', 'seeded', '--seed-decoder pymatching --second-decoder correlated_pymatching')
]

SHOTS = 10000
COMMON_OPTS = f'--mc-max-samples {SHOTS} -p 1e-3 -pp --latency-histogram'

for (folder_prefix, decoder_name, opts) in DECODERS:
    output_folder = f'out/{folder_prefix}_latency'
    os.system(f'mkdir -p {output_folder}')
    for d in [27, 29, 31]:
        output_file = f'{output_folder}/d{d}.out'
        cmd = f'mpirun -np {CORES} ./build/qudec {decoder_name} -d {d} -r {3*d} {COMMON_OPTS} {opts} > {output_file}'
        print(cmd)
        os.system(cmd)
