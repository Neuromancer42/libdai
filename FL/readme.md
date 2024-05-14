# SmartFL with optimized bp

## Setup
get SmartFL from https://github.com/toledosakasa/SMARTFL and set up Defects4j according to its readme

```
cd "Path2SmartFL"
git apply "Path2ThisRepo"/bp.patch
mkdir bptime/ori
mkdir bptime/opt
```

## Generate traces
```
python3 s.py runbp
```

## Test bp time
```
python3 s.py parsebp ori [clean]
python3 s.py parsebp opt [clean]
```

Original bp time is in "Path2SmartFL"/bptime/ori/bp.log  
Optimized bp time is in "Path2SmartFL"/bptime/opt/bp.log  
Add "clean" to remove previous log

## Calculate acceleration rate
```
cd bptime
python3 count.py
```