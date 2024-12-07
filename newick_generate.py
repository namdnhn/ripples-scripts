#!/usr/bin/env python3

import json
import re
from warnings import warn
from collections import defaultdict

with open('asia_data/ncov_gisaid_asia_all-time.json') as f:
    ncov = json.load(f)

# Genes (to which AA mutations are relative):
genePos = {}
for gene,anno in ncov['meta']['genome_annotations'].items():
    if (gene != 'nuc'):
        genePos[gene] = anno['start'] - 1

# Variants and clades

snvRe = re.compile('^([ACGT-])([0-9]+)([ACGT-])$')
snvAaRe = re.compile('^([A-Z*-])([0-9]+)([A-Z*-])$')

newClades = {}
variantCounts = {}
variantAaChanges = {}
samples = []

# Clades from Feb. 2023:
# nColors=23
# head -$nColors ~/github/ncov/defaults/color_schemes.tsv | tail -1 \
# | perl -wne 'chomp; @cols = grep { s/^#//; } split;
#              foreach $hexC (@cols) {
#                if ($hexC =~ m/^(..)(..)(..)$/) {
#                  print sprintf("%d,%d,%d\n", hex($1), hex($2), hex($3));
#                } else { die $hexC; } }'
newCladeColors = { '20H (Beta, V2)':  '94,29,157',
                   '20I (Alpha, V1)': '74,39,178',
                   '20J (Gamma, V3)': '64,58,196',
                   '21A (Delta)' :    '63,82,205',
                   '21I (Delta)':     '65,105,207',
                   '21J (Delta)':     '69,126,203',
                   '21B (Kappa)':     '76,144,192',
                   '21C (Epsilon)':   '85,158,177',
                   '21D (Eta)':       '95,169,160',
                   '21E (Theta)':     '107,177,142',
                   '21F (Iota)':      '121,183,124',
                   '21G (Lambda)':    '137,187,107',
                   '21H (Mu)':        '153,189,93',
                   '21M (Omicron)':   '170,189,82', # B.1.1.529
                   '21K (Omicron)':   '187,188,73', # BA.1
                   '21L (Omicron)':   '203,184,66', # BA.2
                   '22A (Omicron)':   '215,175,62', # BA.4
                   '22B (Omicron)':   '224,162,58', # BA.5
                   '22C (Omicron)':   '230,146,55', # BA.2.12.1
                   '22D (Omicron)':   '230,123,50', # BA.2.75
                   '22E (Omicron)':   '227,97,45',  # BQ.1
                   '22F (Omicron)':   '223,69,40',  # XBB
                   '23A (Omicron)':   '219,40,35',  # XBB.1.5
                   # Grayscale for pre-VoC lineages
                   '19A':       '216,216,216',
                   '19B':       '209,209,209',
                   '20A':       '202,202,202',
                   '20B':       '195,195,195',
                   '20C':       '188,188,188',
                   '20D':       '181,181,181',
                   '20E (EU1)': '174,174,174',
                   '20F':       '167,167,167',
                   '20G':       '160,160,160',
         }

def cladeColorFromName(cladeName, cladeColors):
    color = cladeColors.get(cladeName);
    if (not color):
        color = '0,0,0'
    return color

def cladeFromVariants(name, variants, varStr):
    """Extract bed12 info from an object whose keys are SNV variant names"""
    clade = {}
    snvEnds = []
    # Watch out for back-mutations which make invalid BED because then we have multiple "blocks"
    # at the same position.  Instead, make a back-mutation cancel out the mutation because the
    # mutation is not found at this node.
    changesByPos = defaultdict(list)
    ixsToRemove = []
    for varName in variants:
        m = snvRe.match(varName)
        if (m):
            ref, pos, alt = m.groups()
            prevMut = changesByPos[pos]
            if (prevMut):
                # If multi-allelic, leave the position in the list; but if back-mutation,
                # remove the original mutation.  In either case, don't add this pos again.
                prevIx, prevRef, prevAlt = prevMut
                if (prevAlt == ref and prevRef == alt):
                    ixsToRemove.append(prevIx)
                    changesByPos[pos] = []
            else:
                ix = len(snvEnds)
                changesByPos[pos] = (ix, ref, alt)
                snvEnds.append(int(pos))
        else:
            warn("cladeFromVariants: no match for SNV '%s'" % (varName))
    if ixsToRemove:
        ixsToRemove.sort(reverse=True)
        for ix in ixsToRemove:
            del snvEnds[ix]
    if snvEnds:
        snvEnds.sort()
        snvStarts = [ e-1 for e in snvEnds ]
        snvSizes = [ 1 for e in snvEnds ]
        clade['thickStart'] = min(snvStarts)
        clade['thickEnd'] = max(snvEnds)
        clade['name'] = name
        clade['varSizes'] = snvSizes
        clade['varStarts'] = snvStarts
        clade['varNames'] = varStr
    return clade

def addDatesToClade(clade, numDateAttrs):
    """Add the numeric dates from ncov.json node_attrs.num_date to clade record"""
    if (numDateAttrs):
        clade['dateInferred'] = numDateAttrs['value']
        clade['dateConfMin'] = numDateAttrs['confidence'][0]
        clade['dateConfMax'] = numDateAttrs['confidence'][1]
    else:
        clade['dateInferred'] = clade['dateConfMin'] = clade['dateConfMax'] = ''

def addCountryToClade(clade, countryAttrs):
    """Add country data from ncov.json node_attrs.country to clade"""
    clade['countryInferred'] = countryAttrs['value']
    conf = countryAttrs.get('confidence')
    clade['countryConf'] = ', '.join([ "%s: %0.5f" % (country, conf)
                                       for country, conf in conf.items()]) if conf else ''

def processClade(branch, tag, branchVariants, branchVarStr, clades):
    """If this is the first time we've seen clade, add it to clades"""
    nodeAttrs = branch['node_attrs']
    if (nodeAttrs.get(tag)):
        cladeName = nodeAttrs[tag]['value']
        if (cladeName != 'unassigned' and not cladeName in clades):
            clades[cladeName] = cladeFromVariants(cladeName, branchVariants, branchVarStr)
            addDatesToClade(clades[cladeName], nodeAttrs.get('num_date'))
            if (nodeAttrs.get('country')):
                addCountryToClade(clades[cladeName], nodeAttrs['country'])
            elif (nodeAttrs.get('division')):
                addCountryToClade(clades[cladeName], nodeAttrs['division'])
            clades[cladeName]['topNode'] = branch

def numDateToYmd(numDate):
    """Convert numeric date (decimal year) to integer year, month, day"""
    year = int(numDate)
    isLeapYear = 1 if (year % 4 == 0) else 0
    # Get rid of the year
    numDate -= year
    # Convert to Julian day
    daysInYear = 366 if isLeapYear else 365
    jDay = int(numDate * daysInYear) + 1
    if (jDay > 334 + isLeapYear):
        month, day = 11, (jDay - 334 - isLeapYear)
    elif (jDay > 304 + isLeapYear):
        month, day = 10, (jDay - 304 - isLeapYear)
    elif (jDay > 273 + isLeapYear):
        month, day = 9, (jDay - 273 - isLeapYear)
    elif (jDay > 243 + isLeapYear):
        month, day = 8, (jDay - 243 - isLeapYear)
    elif (jDay > 212 + isLeapYear):
        month, day = 7, (jDay - 212 - isLeapYear)
    elif (jDay > 181 + isLeapYear):
        month, day = 6, (jDay - 181 - isLeapYear)
    elif (jDay > 151 + isLeapYear):
        month, day = 5, (jDay - 151 - isLeapYear)
    elif (jDay > 120 + isLeapYear):
        month, day = 4, (jDay - 120 - isLeapYear)
    elif (jDay > 90 + isLeapYear):
        month, day = 3, (jDay - 90 - isLeapYear)
    elif (jDay > 59 + isLeapYear):
        month, day = 2, (jDay - 59 - isLeapYear)
    elif (jDay > 31):
        month, day = 1, (jDay - 31)
    else:
        month, day = 0, jDay
    return year, month, day

months = ['Jan', 'Feb', 'Mar', 'Apr', 'May', 'Jun', 'Jul', 'Aug', 'Sep', 'Oct', 'Nov', 'Dec'];

def numDateToYmdStr(numDate):
    """Convert decimal year to YY-MM-DD string"""
    if numDate:
        year, month, day = numDateToYmd(numDate)
        return "%02d-%02d-%02d" % (year, month+1, day)
    else:
        return ''

def rUnpackNextstrainTree(branch, parentVariants, parentVarStr):
    """Recursively descend ncov.tree and build data structures for genome browser tracks"""
    # Gather variants specific to this node/branch (if any)
    localVariants = []
    if (branch.get('branch_attrs') and branch['branch_attrs'].get('mutations') and
        branch['branch_attrs']['mutations'].get('nuc')):
        # Nucleotide variants specific to this branch
        for varName in branch['branch_attrs']['mutations']['nuc']:
            if (snvRe.match(varName)):
                localVariants.append(varName)
                if (not variantCounts.get(varName)):
                    variantCounts[varName] = 0;
            else:
                warn("rUnpackNextstrainTree: no snvRe match for '%s'" % (varName))
        # Amino acid variants: figure out which nucleotide variant goes with each
        for geneName in branch['branch_attrs']['mutations'].keys():
            if (geneName != 'nuc'):
                for change in branch['branch_attrs']['mutations'][geneName]:
                    # Get nucleotide coords & figure out which nuc var this aa change corresponds to
                    aaM = snvAaRe.match(change)
                    if (aaM):
                        aaRef, aaPos, aaAlt = aaM.groups()
                        varStartMin = (int(aaPos) - 1) * 3
                        if (genePos.get(geneName)):
                            cdsStart = genePos.get(geneName)
                            varStartMin += cdsStart
                            varStartMax = varStartMin + 2
                            for varName in localVariants:
                                ref, pos, alt = snvRe.match(varName).groups()
                                pos = int(pos) - 1
                                if (pos >= varStartMin and pos <= varStartMax):
                                    variantAaChanges[varName] = geneName + ':' + change
                        else:
                            warn("Can't find start for gene " + geneName)
    # Inherit parent variants
    branchVariants = parentVariants.copy()
    # Add variants specific to this branch (if any)
    for varName in localVariants:
        branchVariants[varName] = 1
    # Make an ordered variant string as David requested: semicolons between nodes,
    # comma-separated within a node:
    branchVarStr = ''
    for varName in localVariants:
        if (len(branchVarStr)):
            branchVarStr += ', '
        branchVarStr += varName
        aaVar = variantAaChanges.get(varName)
        if (aaVar):
            branchVarStr += ' (' + aaVar + ')'
    if (len(parentVarStr) and len(branchVarStr)):
        branchVarStr = parentVarStr + '; ' + branchVarStr
    elif (not len(branchVarStr)):
        branchVarStr = parentVarStr
    processClade(branch, 'clade_membership', branchVariants, branchVarStr, newClades)
    kids = branch.get('children')
    if (kids):
        for child in kids:
            rUnpackNextstrainTree(child, branchVariants, branchVarStr);
    else:
        for varName in branchVariants:
            variantCounts[varName] += 1
        nodeAttrs = branch['node_attrs']
        if (nodeAttrs.get('submitting_lab')):
            lab = nodeAttrs['submitting_lab']['value']
        else:
            lab = ''
        epiNode = nodeAttrs.get('gisaid_epi_isl')
        epiId = epiNode['value'] if epiNode else branch['name']
        numDateNode = nodeAttrs.get('num_date')
        numDate = numDateNode['value'] if numDateNode else ''
        samples.append({ 'id': epiId,
                         'name': branch['name'],
                         'clade': nodeAttrs['clade_membership']['value'],
                         'date': numDateToYmdStr(numDate),
                         'lab': lab,
                         'variants': branchVariants,
                         'varStr': branchVarStr })

rUnpackNextstrainTree(ncov['tree'], {}, '')

def sampleName(sample):
    if sample['id'] != sample['name']:
        return '|'.join([sample['name'], sample['id'], sample['date']])
    else:
        return '|'.join([sample['name'], sample['date']])

sampleCount = len(samples)
sampleNames = [ sampleName(sample)  for sample in samples ]


def sampleIdsFromNode(node, cladeTops=()):
    """Return a list of IDs of all samples found under node."""
    kids = node.get('children')
    if (kids):
        sampleIds = []
        for kid in kids:
            if (kid not in cladeTops):
                sampleIds += sampleIdsFromNode(kid, cladeTops)
    else:
        epiNode = node['node_attrs'].get('gisaid_epi_isl')
        sampleId = epiNode['value'] if epiNode else node['name']
        sampleIds = [sampleId]
    return sampleIds

def sanitizeFileName(filename):
    """Remove or replace characters that cause trouble in filenames"""
    filename = filename.replace('/', '_').replace(' ', '_')
    filename = filename.replace('(', '').replace(')', '')
    filename = filename.replace(',', '');
    return filename

cladeSampleCounts = {}
cladeSampleNames = {}

newCladeTops = [ newClades[cladeName]['topNode'] for cladeName in newClades ]

def rNextstrainToNewick(node, cladeColors, cladeTops=(), parentClade=None, parentVarStr=''):
    """Recursively descend ncov.tree and build Newick tree string of samples to file.
    Exclude nodes in cladeTops."""
    kids = node.get('children')
    if (kids):
        # Make a more concise variant path string than the one we make for the clade track,
        # to embed in internal node labels for Yatish's tree explorations.
        localVariants = []
        if (node.get('branch_attrs') and node['branch_attrs'].get('mutations') and
            node['branch_attrs']['mutations'].get('nuc')):
            # Nucleotide variants specific to this branch
            for varName in node['branch_attrs']['mutations']['nuc']:
                if (snvRe.match(varName)):
                    localVariants.append(varName)
        varStr = '+'.join(localVariants)
        if (len(parentVarStr) and len(varStr)):
            varStr = '$'.join([parentVarStr, varStr])
        elif (not len(varStr)):
            varStr = parentVarStr
        nodeAttrs = node['node_attrs']
        if (nodeAttrs.get('clade_membership')):
            cladeName = nodeAttrs['clade_membership']['value']
        elif (parentClade):
            cladeName = parentClade
        else:
            cladeName = 'unassigned'
        cladeShortened = cladeName.split(' ', 1)[0]
        cladeShortened = cladeShortened.split('/', 1)[0]
        descendants = ','.join([ rNextstrainToNewick(child, cladeColors, cladeTops, cladeName,
                                                     varStr)
                                 for child in kids if child not in cladeTops ])
        treeString = '(' + descendants + ')'
    else:
        nodeAttrs = node['node_attrs']
        epiNode = nodeAttrs.get('gisaid_epi_isl')
        gId = epiNode['value'] if epiNode else node['name']
        name = node['name']
        numDateNode = nodeAttrs.get('num_date')
        date = numDateToYmdStr(numDateNode['value'] if numDateNode else '')
        cladeName = nodeAttrs['clade_membership']['value']
        treeString = sampleName({ 'id': gId, 'name': name, 'date': date })
    return treeString

with open('nextstrain.nh', 'w') as outF:
    outF.write(rNextstrainToNewick(ncov['tree'], newCladeColors) + ';\n')


def newickForClades(clades, cladeColors, cladeTops=()):
    for cladeName in clades:
        filename = 'nextstrain' + sanitizeFileName(cladeName) + '.nh'
        node = clades[cladeName]['topNode']
        with open(filename, 'w') as outF:
            outF.write(rNextstrainToNewick(node, cladeColors, cladeTops) + ';\n')

newickForClades(newClades, newCladeColors, newCladeTops)