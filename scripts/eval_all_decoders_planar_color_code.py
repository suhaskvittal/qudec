import os
cores = 144

PLIST ={
#   'capacity': [0.12, 0.11, 0.1, 0.075, 0.05, 0.025],
#   'pheno': [0.03, 0.025, 0.02, 0.015, 0.01, 0.009],
    'circuit': [0.002, 0.0015, 0.001, 0.00085]
}

EXE = f'mpirun -np {cores} ./build/qudec'
ENSEMBLE_SIZE = 20
SHOTS = 10_000_000
TESSERACT_OPTS = '--tesseract-beam-climbing --tesseract-det-order-bfs --tesseract-num-det-orders 21 --tesseract-det-order-index --tesseract-beam 20'

def run(cmd):
    print(cmd)
    os.system(cmd)

# make folders
for code_type in ['cc488', 'cc666']:
    for f in ['concat_mwpm', 'concat_mwpm_bp_renorm', 'concat_mwpm_ensemble', 'tesseract']:
        for error_model in PLIST:
            os.system(f'mkdir -p out/{f}/{code_type}/{error_model}')

for (error_model, plist) in PLIST.items():
    # resolve model option
    model_opt = ''
    if error_model == 'capacity':
        model_opt = '-cap'
    elif error_model == 'pheno':
        model_opt = '-pheno'

    # resolve flag protocol option
    flag_opt_1, flag_opt_2 = '', ''
    if error_model == 'circuit':
        flag_opt_1 = '--concat-flag-renormalization'
        flag_opt_2 = '--concat-bp'

    for p in plist:
        for code_type in ['cc488', 'cc666']:
            for d in [5, 7, 9, 11, 13]:
                folder_suffix, filename = f'{code_type}/{error_model}', f'd{d}_p{p:.2e}.out'

                concat_mwpm_flag_opt = '--concat-flag-renormalization' if error_model == 'circuit' else ''
                ensemble_flag_opt = '--concat-bp' if error_model == 'circuit' else ''

                ## Concat-MWPM
                run(f'{EXE} concat_mwpm {d} -p {p} --seed 0 --mc-max-samples {SHOTS} --code {code_type} {model_opt} {concat_mwpm_flag_opt} -pp > out/concat_mwpm/{folder_suffix}/{filename}')
                ## Concat-MWPM with ensembling
                run(f'{EXE} concat_mwpm {d} -p {p} --seed 0 --mc-max-samples {SHOTS} --code {code_type} {model_opt} --concat-ensemble {ENSEMBLE_SIZE} {ensemble_flag_opt} -pp > out/concat_mwpm_ensemble/{folder_suffix}/{filename}')

                ## Concat-MWPM with BP renormalization
                if error_model == 'circuit':
                    run(f'{EXE} concat_mwpm {d} -p {p} --seed 0 --mc-max-samples {SHOTS} --code {code_type} {model_opt} --concat-bp -pp > out/concat_mwpm_bp_renorm/{folder_suffix}/{filename}')

                ## Tesseract
                if d < 11:
                    run(f'{EXE} tesseract {d} -p {p} --seed 0 --mc-max-samples {SHOTS} --code {code_type} {model_opt} {TESSERACT_OPTS} -pp > out/tesseract/{folder_suffix}/{filename}')

