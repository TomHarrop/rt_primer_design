# -*- coding: utf-8 -*-
#
############################
# Design real-time primers #
############################

# Imports

import os
import argparse
from primerDesignFunctions import runBlast
from subprocess import call
import csv
import urllib

# command-line options
parser = argparse.ArgumentParser(
    description='Design primers for real-time PCR from a list of LOC IDs.')
parser.add_argument('--output', '-o', help='Output directory', type=str,
                    required=True, dest='output')
parser.add_argument('--input', '-i', help='List of LOC IDs (one per line)',
                    type=str, required=True, dest='primersList')
args = parser.parse_args()
outdir = args.output
primersList = args.primersList

#############
# variables #
#############

# primer-BLAST parameters.
parameters = {
    'OVERLAP_5END': '7',
    'OVERLAP_3END': '4',
    'PRIMER_PRODUCT_MIN': '75',
    'PRIMER_PRODUCT_MAX': '180',
    'PRIMER_NUM_RETURN': '10',
    'PRIMER_MIN_TM': '55',
    'PRIMER_OPT_TM': '60',
    'PRIMER_MAX_TM': '65',
    'PRIMER_MAX_DIFF_TM': '5',
    'PRIMER_ON_SPLICE_SITE': '0',
    'SEARCHMODE': '0',
    'SPLICE_SITE_OVERLAP_5END': '7',
    'SPLICE_SITE_OVERLAP_3END': '4',
    'SPAN_INTRON': 'checked',
    'MIN_INTRON_SIZE': '50',
    'MAX_INTRON_SIZE': '1000000',
    'SEARCH_SPECIFIC_PRIMER': 'on',
    'EXCLUDE_ENV': 'off',
    'EXCLUDE_XM': 'off',
    'TH_OLOGO_ALIGNMENT': 'off',
    'TH_TEMPLATE_ALIGNMENT': 'off',
    'ORGANISM': 'Oryza sativa Japonica Group (taxid:39947)',
    'PRIMER_SPECIFICITY_DATABASE': 'refseq_rna',
    'TOTAL_PRIMER_SPECIFICITY_MISMATCH': '1',
    'PRIMER_3END_SPECIFICITY_MISMATCH': '1',
    'MISMATCH_REGION_LENGTH': '5',
    'TOTAL_MISMATCH_IGNORE': '7',
    'PRODUCT_SIZE_DEVIATION': '200',
    'ALLOW_TRANSCRIPT_VARIANTS': 'on',
    'HITSIZE': '50000',
    'EVALUE': '30000',
    'WORD_SIZE': '7',
    'MAX_CANDIDATE_PRIMER': '1000',
    'PRIMER_MIN_SIZE': '15',
    'PRIMER_OPT_SIZE': '20',
    'PRIMER_MAX_SIZE': '25',
    'PRIMER_MIN_GC': '45',
    'PRIMER_MAX_GC': '55',
    'GC_CLAMP': '2',
    'NUM_TARGETS_WITH_PRIMERS': '1000',
    'NUM_TARGETS': '20',
    'MAX_TARGET_PER_TEMPLATE': '100',
    'POLYX': '3',
    'SELF_ANY': '3',
    'SELF_END': '1',
    'PRIMER_MAX_END_STABILITY': '9',
    'PRIMER_MAX_END_GC': '5',
    'PRIMER_MAX_TEMPLATE_MISPRIMING_TH': '40.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING_TH': '70.00',
    'PRIMER_MAX_SELF_ANY_TH': '45.0',
    'PRIMER_MAX_SELF_END_TH': '35.0',
    'PRIMER_PAIR_MAX_COMPL_ANY_TH': '45.0',
    'PRIMER_PAIR_MAX_COMPL_END_TH': '35.0',
    'PRIMER_MAX_HAIRPIN_TH': '24.0',
    'PRIMER_MAX_TEMPLATE_MISPRIMING': '12.00',
    'PRIMER_PAIR_MAX_TEMPLATE_MISPRIMING': '24.00',
    'PRIMER_PAIR_MAX_COMPL_ANY': '8.00',
    'PRIMER_PAIR_MAX_COMPL_END': '3.00',
    'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
    'NO_SNP': 'off',
    'LOW_COMPLEXITY_FILTER': 'on',
    'MONO_CATIONS': '50.0',
    'DIVA_CATIONS': '1.5',
    'CON_ANEAL_OLIGO': '50.0',
    'CON_DNTPS': '0.6',
    'SALT_FORMULAR': '1',
    'TM_METHOD': '1',
    'PRIMER_INTERNAL_OLIGO_MIN_SIZE': '18',
    'PRIMER_INTERNAL_OLIGO_OPT_SIZE': '20',
    'PRIMER_INTERNAL_OLIGO_MAX_SIZE': '27',
    'PRIMER_INTERNAL_OLIGO_MIN_TM': '57.0',
    'PRIMER_INTERNAL_OLIGO_OPT_TM': '60.0',
    'PRIMER_INTERNAL_OLIGO_MAX_TM': '63.0',
    'PRIMER_INTERNAL_OLIGO_MAX_GC': '80.0',
    'PRIMER_INTERNAL_OLIGO_OPT_GC_PERCENT': '50',
    'PRIMER_INTERNAL_OLIGO_MIN_GC': '20.0',
    'PICK_HYB_PROBE': 'off',
    'NEWWIN': 'off',
    'SHOW_SVIEWER': 'false',
}

########
# CODE #
########

# call R script to get refseq ids
call(['bin/getRefSeqIDs.R', primersList, outdir])

print(('\nRefSeq IDs for genes in excludedRecords.csv may be manually added'
       '\nto RAP.MSU.refseq.csv to process them.\n\nPress enter to continue'
       '\nprocessing RAP.MSU.refseq.csv'))
input()

# 2. Read output from getRefSeqIDs. Make a dict of MSU ID : RefSeq ID

RapMsuRefSeq = csv.reader(open(outdir + '/RAP.MSU.refseq.csv'))
next(RapMsuRefSeq, None)
MsuRefSeq = {}
for row in RapMsuRefSeq:
    MsuRefSeq[row[0]] = row[1]

# 3. initial call to runBlast function

ResultsPages = runBlast(MsuRefSeq, parameters=parameters)

# primers that worked first time
strictPrimers = {}
for key in ResultsPages:
    if (not ResultsPages[key].exceptions and
            not ResultsPages[key].offTargets and
            not ResultsPages[key].noPrimersFound):
        ResultsPages[key].finalStatus = 'OK'
        strictPrimers[key] = ResultsPages[key]

for key in strictPrimers:
    ResultsPages.pop(key)

# 4. iterate the BLAST search:

# a. first, remove genes without introns
intronlessReruns = {}
intronlessPrimers = {}
for key in ResultsPages:
    if ResultsPages[key].exceptions:
        intronlessReruns[key] = ResultsPages[key].RefSeq
        intronlessPrimers[key] = ResultsPages[key]
        intronlessPrimers[key].finalStatus = 'failed_no_introns'
for key in intronlessReruns:
    ResultsPages.pop(key)
if len(intronlessReruns) > 0:
    print('Removed ' + str(len(intronlessReruns)) +
          ' genes without introns.\n')

# =============================================================================
# # For now rerunning the intronless genes is not possible because the BLAST
# # server is ignoring the SPAN_INTRON input. Later, output URLs for the
# # intronless genes to allow easy manual designing.
#
# intronlessPrimers = runBlast(intronlessReruns, intronlessParams)
# intronlessStrict = {}
# for key in intronlessPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         intronlessPrimers[key].finalStatus = 'INTRONLESS_strict'
#         intronlessStrict[key] = intronlessPrimers[key]
# for key in intronlessStrict:
#     intronlessReruns.pop(key)
#
# intronlessParams['GC_CLAMP'] = '1'
# intronlessPrimers = runBlast(intronlessReruns, intronlessParams)
# intronlessGC1 = {}
# for key in intronlessPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         intronlessPrimers[key].finalStatus = 'INTRONLESS_GC1'
#         intronlessGC1[key] = intronlessPrimers[key]
# for key in intronlessGC1:
#     intronlessReruns.pop(key)
#
# intronlessParams['PRIMER_MIN_GC'] = '40'
# intronlessParams['PRIMER_MAX_GC'] = '60'
# intronlessPrimers = runBlast(intronlessReruns, intronlessParams)
# intronlessGCDev = {}
# for key in intronlessPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         intronlessPrimers[key].finalStatus = 'INTRONLESS_GCDev'
#         intronlessGCDev[key] = intronlessPrimers[key]
# for key in intronlessGCDev:
#     intronlessReruns.pop(key)
#
# intronlessParams['PRIMER_MIN_TM'] = '52'
# intronlessPrimers = runBlast(intronlessReruns, intronlessParams)
# intronlessLowTM = {}
# for key in intronlessPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         intronlessPrimers[key].finalStatus = 'INTRONLESS_GCDev'
#         intronlessLowTM[key] = intronlessPrimers[key]
# for key in intronlessLowTM:
#     intronlessReruns.pop(key)
#
# intronlessParams['SELF_ANY'] = '5'
# intronlessParams['SELF_END'] = '2'
# intronlessPrimers = runBlast(intronlessReruns, intronlessParams)
# intronlessDimers = {}
# for key in intronlessPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         intronlessPrimers[key].finalStatus = 'INTRONLESS_GCDev'
#         intronlessDimers[key] = intronlessPrimers[key]
# for key in intronlessDimers:
#     intronlessReruns.pop(key)
# =============================================================================

# b. relax perameters and re-iterate
relaxedParams = parameters.copy()
relaxedReruns = {}
for key in ResultsPages:
    relaxedReruns[key] = ResultsPages[key].RefSeq

# i. GC_CLAMP
print('Found ' + str(len(strictPrimers)) + ' strict primer pairs.\n')
print('Trying again with relaxed parameters: GC Clamp\n')

relaxedParams['GC_CLAMP'] = '1'
if relaxedReruns:
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedGC1 = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'GC1'
        relaxedGC1[key] = relaxedPrimers[key]
for key in relaxedGC1:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)
relaxedParams['GC_CLAMP'] = '0'
if relaxedReruns:
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedGC0 = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'GC0'
        relaxedGC0[key] = relaxedPrimers[key]
for key in relaxedGC0:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)

# ii. Primer GC content
if len(relaxedGC1) + len(relaxedGC0) > 0:
    print('Found another ' + str(len(relaxedGC1) + len(relaxedGC0)) +
          '  primer pairs.\n')

relaxedParams['PRIMER_MIN_GC'] = '35'
relaxedParams['PRIMER_MAX_GC'] = '65'
if relaxedReruns:
    print('Trying again with relaxed parameters: GC Content\n')
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedGCCont = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'GC_Content'
        relaxedGCCont[key] = relaxedPrimers[key]

for key in relaxedGCCont:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)

# iii. Primer TM
if relaxedGCCont:
    print('Found another ' + str(len(relaxedGCCont)) + '  primer pairs.\n')

relaxedParams['PRIMER_MIN_TM'] = '52'
if relaxedReruns:
    print('Trying again with relaxed parameters: Primer TM\n')
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedTM = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'Low_TM'
        relaxedTM[key] = relaxedPrimers[key]
for key in relaxedTM:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)

# iv. Primer complementarity
if relaxedTM:
    print('Found another ' + str(len(relaxedTM)) + '  primer pairs.\n')

relaxedParams['SELF_ANY'] = '5'
relaxedParams['SELF_END'] = '2'
if relaxedReruns:
    print('Trying again with relaxed parameters: Self complementarity\n')
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedDimers = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'Potential_Dimers'
        relaxedDimers[key] = relaxedPrimers[key]
for key in relaxedDimers:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)

# v. Near-defaults (use with caution)
if relaxedDimers:
    print('Found another ' + str(len(relaxedDimers)) + '  primer pairs.\n')

relaxedParams['SELF_ANY'] = '8'
relaxedParams['SELF_END'] = '3'
if relaxedReruns:
    print('Trying again with DEFAULT parameters: Self complementarity\n'
          'n.b. These primers may not work well.\n')
    relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
else:
    relaxedPrimers = {}

relaxedBad = {}
for key in relaxedPrimers:
    if (not relaxedPrimers[key].exceptions and
            not relaxedPrimers[key].offTargets and
            not relaxedPrimers[key].noPrimersFound):
        relaxedPrimers[key].finalStatus = 'Probable_Dimers'
        relaxedBad[key] = relaxedPrimers[key]
for key in relaxedBad:
    relaxedReruns.pop(key)
    ResultsPages.pop(key)

# =============================================================================
# # vi. Without repeat filtering (doesn't work)
# relaxedParams['LOW_COMPLEXITY_FILTER'] = 'off'
# relaxedPrimers = runBlast(relaxedReruns, relaxedParams)
# relaxedNoFilter = {}
# for key in relaxedPrimers:
#    if (not relaxedPrimers[key].exceptions and
#            not relaxedPrimers[key].offTargets and
#            not relaxedPrimers[key].noPrimersFound):
#         relaxedPrimers[key].finalStatus = 'No_repeat_filter'
#         relaxedNoFilter[key] = relaxedPrimers[key]
# for key in relaxedNoFilter:
#     relaxedReruns.pop(key)
#     ResultsPages.pop(key)
# =============================================================================

for key in ResultsPages:
    if ResultsPages[key].exceptions:
        ResultsPages[key].finalStatus = 'unknown_error'
    if ResultsPages[key].offTargets:
        ResultsPages[key].finalStatus = 'no_specific_primers'
    if ResultsPages[key].noPrimersFound:
        ResultsPages[key].finalStatus = 'primer_quality_too_low'

primerSets = [strictPrimers, relaxedGC1, relaxedGC0, relaxedBad, relaxedGCCont,
              relaxedDimers]
primerStrings = ['strictPrimers', 'relaxedGC1', 'relaxedGC0', 'relaxedBad',
                 'relaxedGCCont', 'relaxedDimers']
failedSets = [ResultsPages, intronlessPrimers]
failedStrings = ['noPrimersFound', 'noIntronsFound']

primerNo = 0
for primerSet in primerSets:
    primerNo += len(primerSet)

print('\nFound ' + str(primerNo) + ' primer pair(s).\nWriting reports...')

with open(outdir + '/primerSummary.csv', 'w') as file:
    file.write('MSU.ID,RefSeqID,Status,PrimerF,PrimerF.TM,PrimerR,PrimerR.TM,'
               'ProductSize,IntronSize\n')
    for primerSet in primerSets:
        for key in primerSet:
            primerSet[key].parsePrimers()
            primerSet[key].replaceCssLinks()
            file.write(primerSet[key].csvLine() + '\n')
    for primerSet in failedSets:
        for key in primerSet:
            file.write('{0},{1},{2},\n'.format(primerSet[key].LOC,
                       primerSet[key].RefSeq, primerSet[key].finalStatus))


# 5. print files to sub-folders of working directory
for i in range(len(primerStrings)):
    if primerSets[i]:
        subdir = outdir + "/" + primerStrings[i]
        if not os.path.isdir(subdir):
            os.mkdir(path=subdir)
        for key in primerSets[i]:
            primerSets[i][key].printFile(subdir=subdir)

# 6. print URLs for manual chase-up into subfolder of working directory
# need to work out how to format the URL without posting it.
intronlessParams = parameters.copy()
intronlessParams['SPAN_INTRON'] = 'off'
for i in range(len(failedStrings)):
    if failedSets[i]:
        subdir = outdir + "/FAILED" + failedStrings[i]
        if not os.path.isdir(subdir):
            os.mkdir(path=subdir)

        with open(subdir + '/links.html', 'w') as file:
            file.write('''<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.01//EN"'''
                       '''\n    "http://www.w3.org/TR/html4/strict.dtd">\n'''
                       '''<html lang="en">\n  '''
                       '''<head>\n    '''
                       '''<meta http-equiv="content-type"'''
                       ''' content="text/html; charset=utf-8">\n'''
                       '''    <title>Links for primer design</title>\n'''
                       '''    <link rel="stylesheet" type="text/css"''
                       '''' href="style.css">\n    '''
                       '''<script type="text/javascript" src="script.js">'''
                       '''</script>\n  </head>\n  <body>\n<p>\n''')
            for key in failedSets[i]:
                pasteParams = intronlessParams.copy()
                pasteParams['INPUT_SEQUENCE'] = failedSets[i][key].RefSeq
                file.write(
                    '''<a href="http://www.ncbi.nlm.nih.gov/tools/'''
                    '''primer-blast/index.cgi?''' +
                    urllib.parse.urlencode(pasteParams) + '''">''' +
                    failedSets[i][key].LOC + '''</a><br />\n''')
            file.write('\n</p>\n</body>\n</html>')

print('\nDone.\n')
