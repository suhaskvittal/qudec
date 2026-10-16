import os
cores = 64
for p in [7e-4, 8.5e-4, 1e-3, 1.5e-3, 2e-3]:
    # pymatching:
    os.system(f'python3 scripts/eval_decoder_uniform.py pymatching {p} {cores} auto')
    # correlated pm:
    os.system(f'python3 scripts/eval_decoder_uniform.py pymatching {p} {cores} correlated_pymatching --correlated')
    # uf:
    os.system(f'python3 scripts/eval_decoder_uniform.py uf {p} {cores} auto')
    # seeded
    os.system(f'python3 scripts/eval_decoder_uniform.py seeded {p} {cores} auto --second-decoder libra --seed-gap-threshold-db 40')
