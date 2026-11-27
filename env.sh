# sourced by run_*.sh: single thread, fixed hash seed, low priority
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONHASHSEED=0
PY=${PY:-python3}
