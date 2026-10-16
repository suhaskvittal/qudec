from sys import argv
import os

decoder = argv[1]
p = float(argv[2])
cores = argv[3]
output_folder = argv[4]
if len(argv) > 5:
    cli_args = ' '.join(argv[5:])
else:
    cli_args = ''

seed = 0
if output_folder == 'auto':
    output_folder = f'out/{decoder}_sc_error_rates_p{p:.2e}'
else:
    output_folder = f'out/{output_folder}_sc_error_rates_p{p:.2e}'

# (d, shots)
comb = [
    (3, 1_000_000),
    (5, 1_000_000),
    (7, 10_000_000),
    (9, 50_000_000),
    (11, 100_000_000),
    (13, 500_000_000),
]


os.system(f'mkdir -p {output_folder}')
for (d, shots) in comb:
    os.system(f'mkdir -p {output_folder}/d{d}')
    # write cli and other info to setup.log
    with open(f'{output_folder}/d{d}/setup.log', 'w') as wr:
        wr.write(f'''p = {p}
d = {d}
shots = {shots}
seed = {seed}
cores = {cores}
cli_args = {cli_args}
''')
    cmd = f'mpirun -np {cores} ./build/qudec {decoder} {d} -p {p} --seed {seed} --mc-max-samples {shots} {cli_args} > {output_folder}/d{d}/output.log'
    print(cmd)
    os.system(cmd)
