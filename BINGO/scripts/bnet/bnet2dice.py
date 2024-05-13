#!/usr/bin/env python3

# Given a bnet.out file produced by cons_all2bnet.py, and a ruleProb.txt file mapping each rule to its probability of
# firing, this script produces a factorGraph.fg file accepted by LibDAI.
# ./bnet2fg.py ruleProb.txt 0.99 < named_bnet.out > factorGraph.fg 2> bnet2fg.log

import logging
import math
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

outLines = [ """network factor
{
}"""]
logging.info(outLines[0])
logging.info('')
print(outLines[0])
print('')

factorLines = bnetLines[1:(numVars + 1)]
assert len(factorLines) == numVars
for varIndex in range(numVars):
    outLines = [f"# Print node {varIndex}"]
    outLines.append(f"variable \"n{varIndex}\"")
    outLines.append("{")
    outLines.append( '    type discrete[2] { \"false\" \"true\" };')
    outLines.append("}")
    # Print.
    for line in outLines:
        line = str(line)
        logging.info(line)
        if not line.startswith('#'): print(line)

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

    outLines.append('# Factor {0} of {1}. Finished printing {2}% of factors.'.format(varIndex, numVars, 100 * varIndex / numVars))
    outLines.append('# {0}'.format(line))

    # The first line contains the clause
    outLines.append(f"probability ( \"n{varIndex}\" | " + "".join([f"\"n{i}\" " for i in parents]) + " )")
    # start the bracket
    outLines.append('{')
    if factorType == '*':
        # The fourth line contains the number of nonzero entries in the factor table.
        tableSize = int(math.pow(2, numParents))

        # The rest of the lines contain these nonzero entries;
        # each line consists of a table index, followed by the value corresponding to that table index.
        # The most difficult part is getting the indexing right.
        # The convention that is used is that the left-most variables
        # cycle through their values the fastest
        # (similar to MatLab indexing of multidimensional arrays).
        p = ruleProbs[ruleName] if (ruleName != None and ruleName in ruleProbs) else \
                      1.0 if ruleName == 'Rnarrow' else \
                      defaultProbability
        q = 0.0
        outLines.append(f"    default {1-q} {q};")
        outLines.append( "    (" + ", ".join(["\"true\"" for _ in parents]) + f") {1-p} {p};")
    else:
        # The fourth line contains the number of nonzero entries in the factor table.
        tableSize = int(math.pow(2, numParents))
        
        p = 1.0
        q = 0.0

        outLines.append(f"    default {q} {1-q};")
        outLines.append( "    (" + ", ".join(["\"false\"" for _ in parents]) + f") {p} {1-p};")

    # close the bracket.
    outLines.append('}')

    # Print.
    for line in outLines:
        line = str(line)
        logging.info(line)
        if not line.startswith('#'): print(line)
