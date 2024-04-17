#!/usr/bin/env python3

# Given a bnet.out file produced by cons_all2bnet.py, and a ruleProb.txt file mapping each rule to its probability of
# firing, this script produces a factorGraph.fg file accepted by LibDAI.
# ./bnet2fg.py ruleProb.txt 0.99 < named_bnet.out > factorGraph.fg 2> bnet2fg.log

import logging
import math
import subprocess
import sys

ruleProbFileName = sys.argv[1]
defaultProbability = float(sys.argv[2])

logging.basicConfig(level=logging.INFO, \
                    format="[%(asctime)s] %(levelname)s [%(name)s.%(funcName)s:%(lineno)d] %(message)s", \
                    datefmt="%H:%M:%S")

########################################################################################################################
# 1. Accept input

# Read bayesian network
bnetLines = [ line.strip() for line in sys.stdin ]

# Load rule probabilities
ruleProbs = [ line.strip().split(': ') for line in open(ruleProbFileName) ]
ruleProbs = { line[0]: float(line[1]) for line in ruleProbs }

########################################################################################################################
# 2. Compute output

# https://staff.fnwi.uva.nl/j.m.mooij/libDAI/doc/fileformats.html
# It starts with a line containing the number of factors in that graph, followed by an empty line.
numVars = int(bnetLines[0])

outLines = [ numVars, '' ]
logging.info(numVars)
logging.info('')
print(numVars)
print('')

factorLines = bnetLines[1:(numVars + 1)]
assert len(factorLines) == numVars
for varIndex in range(numVars):
    outLines = []
    # Then all factors are specified, using one block for each factor,

    line = factorLines[varIndex]
    components = [ c.strip() for c in line.split(' ') ]

    factorType = components[0] # '*' or '+'
    assert factorType == '*' or factorType == '+'
    ruleName = components[1] if factorType == '*' else None
    numParents = int(components[2]) if factorType == '*' else int(components[1])
    parents = [ int(p) for p in components[3:] ] if factorType == '*' else [ int(p) for p in components[2:] ]

    # Each causal factor starts with a head variable
    outLines.append(varIndex)
    # The second line contains the factor type and the probability
    if factorType == '*':
        probability = ruleProbs[ruleName] if (ruleName != None and ruleName in ruleProbs) else \
                      1.0 if ruleName == 'Rnarrow' else \
                      defaultProbability
        outLines.append('*'+str(probability))
    else:
        outLines.append('+')
    # The third line is the number of parents
    outLines.append(numParents)
    # The forth line are the parents
    outLines.append(' '.join([str(p) for p in parents]))

    # where the factors are seperated by empty lines.
    outLines.append('')

    # Print.
    for line in outLines:
        line = str(line)
        logging.info(line)
        if not line.startswith('#'): print(line)
