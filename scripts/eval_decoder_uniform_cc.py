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
    output_folder = f'out/{decoder}_cc_capacity_error_rates_p{p:.2e}'
else:
    output_folder = f'out/{output_folder}_cc_capacity_error_rates_p{p:.2e}'

# (d, shots)
SHOTS = 10_000_000
comb = [
    (d, SHOTS) for d in range(5, 31+1, 2)
]
ETH = 5 / SHOTS


os.system(f'mkdir -p {output_folder}')
previous_ler = None
previous_distance = None
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
    if previous_ler is not None and previous_ler < ETH:
        reason = f'previous evaluated distance d{previous_distance} had zero logical error rate'
        with open(f'{output_folder}/d{d}/output.log', 'w') as wr:
            wr.write(f'''Logical error rate: 0
Dummy output: run skipped because {reason}.
Timing data: invalid (no decoding was performed).
''')
        print(f'Skipping d{d}: {reason}')
        continue

    cmd = f'mpirun -np {cores} ./build/qudec {decoder} {d} -p {p} --seed {seed} --mc-max-samples {shots} {cli_args} --code-capacity > {output_folder}/d{d}/output.log'
    print(cmd)
    status = os.system(cmd)
    previous_ler = None
    previous_distance = d
    # Only a successful run with an unambiguous rate can trigger skipping.
    if status == 0:
        with open(f'{output_folder}/d{d}/output.log') as rd:
            rates = [line.split(':', 1)[1].strip() for line in rd
                     if line.startswith('Logical error rate:')]
        if len(rates) == 1:
            try:
                previous_ler = float(rates[0])
            except ValueError:
                pass
