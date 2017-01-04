#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import rtPrimerDesign.functions as functions


#############
# PARSE CLI #
#############

# FIXME: should allow the organism to be passed in

# FIXME: expect a list/dict of RefSeq IDs from cli

#     LOC number for testing : LOC_Os05g41760
#         RAP-ID for testing : Os05g0497200
# RefSeq mRNA ID for testing : NM_001062476
test_ref_seq = 'NM_001062476'

##############
# PARAMETERS #
##############

strict_parameters = {
    'PRIMER_PRODUCT_MIN': '70',
    'PRIMER_PRODUCT_MAX': '180',
    'PRIMER_NUM_RETURN': '10',
    'PRIMER_MIN_TM': '55.0',
    'PRIMER_OPT_TM': '60.0',
    'PRIMER_MAX_TM': '65.0',
    'PRIMER_MAX_DIFF_TM': '5',
    'MIN_INTRON_SIZE': '0',
    'MAX_INTRON_SIZE': '1000000',
    'PRIMER_SPECIFICITY_DATABASE': 'refseq_mrna',
    'EXCLUDE_ENV': 'on',
    'ORGANISM': 'Oryza sativa Japonica Group (taxid:39947)',
    'TOTAL_MISMATCH_IGNORE': '7',
    'ALLOW_TRANSCRIPT_VARIANTS': 'on',
    'MAX_CANDIDATE_PRIMER': '1000',
    'PRIMER_MIN_GC': '45.0',
    'PRIMER_MAX_GC': '55.0',
    'GC_CLAMP': '2',
    'POLYX': '3',
    'SELF_ANY': '3.00',
    'SELF_END': '1.00',
    'SEARCH_SPECIFIC_PRIMER': 'on',
    'SHOW_SVIEWER': 'on',
    'UNGAPPED_BLAST': 'on',
    'LOW_COMPLEXITY_FILTER': 'on',
    'SHOW_SVIEWER': 'on',
    'SPAN_INTRON': 'on'
}

#############
# FUNCTIONS #
#############

# define the looping syntax here.

# submit a BLAST job for the goi
blast_responses = {}
blast_responses[test_ref_seq] = functions.primerBlastResults(
    test_ref_seq, parameters)

# wait for job to finish
blast_responses[test_ref_seq].checkRunning()
while blast_responses[test_ref_seq].running:
    # wait here
    blast_responses[test_ref_seq].pollResults()
    blast_responses[test_ref_seq].checkRunning()

# check for similar targets
blast_responses[test_ref_seq].user_seqloc
blast_responses[test_ref_seq].check_similar_templates()

# TRY_USER_GUIDE=yes&USER_SEQLOC=ref|XM_015782903.1|?455?1369&FIRST_USER_RID=6BH3N8SH014
# re-blast in user guided mode
post_parameters = blast_responses[test_ref_seq].blast_parameters.copy()
post_parameters['TRY_USER_GUIDE'] = 'yes'
post_parameters['USER_SEQLOC'] = blast_responses[test_ref_seq].user_seqloc
blast_response = requests.get(
                blast_responses[test_ref_seq].blastUrl,
                params=post_parameters)
BeautifulSoup(blast_response.content, 'lxml')



