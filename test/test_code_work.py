#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import rtPrimerDesign.functions as functions
import time

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

strict_parameters_intronless = strict_parameters.copy()
strict_parameters_intronless.pop('SPAN_INTRON')


# function to run primer blast and wait for it to finish
def run_primer_blast(refseq_id, blast_parameters):
    # run a BLAST search
    blast_result = functions.primerBlastResults(
        RefSeq=refseq_id,
        blast_parameters=blast_parameters)
    # wait for job to finish
    blast_result.checkRunning()
    while blast_result.running:
        time.sleep(10)
        blast_result.pollResults()
        blast_result.checkRunning()
    # check for exon/exon junction
    blast_result.check_introns()
    # check for similar sequences and re-run if required
    blast_result.check_similar_templates()
    blast_result.checkRunning()
    while blast_result.running:
        time.sleep(10)
        blast_result.pollResults()
        blast_result.checkRunning()
    # return primerBlastResults
    return(blast_result)

# gene of interest
# LOC_Os02g41460,XM_015768655,OS02G0623400,G1L3
# LOC_Os05g41760 Os05g0497200, XM_015782903.1, NM_001062476
test_refseq = 'XM_015782903'
test_blast_result = run_primer_blast(test_refseq, strict_parameters)

# re-run without intron span
if test_blast_result.no_intron:
    test_blast_result = run_primer_blast(
        test_refseq,
        strict_parameters_intronless)

# re-run with lower stringency

