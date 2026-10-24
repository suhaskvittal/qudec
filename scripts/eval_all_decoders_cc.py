import os
cores = 64

PLIST = [
    0.1,
    0.075,
    0.05,
    0.025,
    0.01
]

for p in PLIST:
    for code_type in ['cc488', 'cc666']:
    #   # restriction
    #   os.system(f'python3 scripts/eval_decoder_uniform_cc.py restriction {p} {cores} auto')

        # concat mwpm
#       os.system(f'python3 scripts/eval_decoder_uniform_cc.py concat_mwpm {p} {cores} concat_mwpm_{code_type} --code {code_type}')
        # concat mwpm with complementary matching
#       os.system(f'python3 scripts/eval_decoder_uniform_cc.py concat_mwpm {p} {cores} concat_mwpm_compl_{code_type} --concat-boundary-gap --code {code_type}')
        # concat mwpm with ensembling
        os.system(f'python3 scripts/eval_decoder_uniform_cc.py concat_mwpm {p} {cores} concat_mwpm_ensembling_{code_type} --concat-ensemble 20 --code {code_type}')

