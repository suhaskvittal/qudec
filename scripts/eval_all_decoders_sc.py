import os
cores = 32

PLIST = [
    8e-4,
    9e-4,
    1e-3,
    1.25e-3,
    1.5e-3,
    1.75e-3,
    2e-3,
    2.5e-3,
    3e-3
]

for p in PLIST:
    # pymatching:
    os.system(f'python3 scripts/eval_decoder_sc.py pymatching {p} {cores} auto')
    # correlated pm:
    os.system(f'python3 scripts/eval_decoder_sc.py pymatching {p} {cores} correlated_pymatching --correlated')
    # uf:
    os.system(f'python3 scripts/eval_decoder_sc.py uf {p} {cores} auto')
    # seeded
    os.system(f'python3 scripts/eval_decoder_sc.py seeded {p} {cores} auto --second-decoder libra --seed-gap-threshold-db 40')
