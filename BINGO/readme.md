# Bingo with optimized bp

## Setup

1. Fetch pre-generated benchmark files and extracted to somewhere

```
tar xfv pjbench.tar.gz -C <path-to-pjbench>
```

2. Set basic environment variable to the location of libdai

```
export LIBDAI=<path-to-libdai>
```

3. Build the bingo wrapper and do other necessary settings

```
source $LIBDAI/BINGO/init.sh <path-to-pjbench>
```

## Run testcases

```
run-opt-bingo.sh <PROJ> ...
```

`<PROJ>` can be one or several of the benchmarks: `hedc | ftp | weblech| jspider | avrora | luindex | sunflow | xalan`

The results would be in the directory `$BINGO_RESULTS`. You can set this environment before run `init.sh`, or it will be `/tmp/results` by default.

The data used in the evaluation is listed in `$BINGO_RESULTS/datarace/<PROJ>/output/exact/stats.<causal|orig>.txt`

*If you want to re-run the test cases, you have to clean up the medium files by `bash $LIBDAI/BINGO/cleanup.sh` (Note: output results would be deleted too!!!)*
