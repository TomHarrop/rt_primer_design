#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################
# DEPENDENCIES #
################

import rtPrimerDesign.functions as functions
import time
import joblib
import tompytools
import sys


#############
# FUNCTIONS #
#############

# function to run primer blast and wait for it to finish
def run_primer_blast(
        refseq_id,
        blast_parameters,
        status,
        wait_seconds=60,
        verbose=False):
    if verbose:
        tompytools.generate_message('Running BLAST for %s' % refseq_id)
        print(
            'RefSeq=%s\t\t\t'
            'blast_parameters=%s\t\t\t'
            'status=%s' %
            (refseq_id, blast_parameters, status))
    # run a BLAST search
    blast_result = functions.primerBlastResults(
        RefSeq=refseq_id,
        blast_parameters=blast_parameters,
        status=status)
    # wait for job to finish
    while blast_result.running:
        if verbose:
            tompytools.generate_message(
                '%s: Waiting %i seconds for BLAST' % (refseq_id, wait_seconds))
        time.sleep(wait_seconds)
        if verbose:
            tompytools.generate_message('%s.pollResults()' % refseq_id)
        blast_result.pollResults()
    # check for exon/exon junction
    if verbose:
        tompytools.generate_message('%s.check_introns()' % refseq_id)
    blast_result.check_introns()
    # check for similar sequences and re-run if required
    if verbose:
        tompytools.generate_message('%s.check_similar_templates()' % refseq_id)
    blast_result.check_similar_templates()
    while blast_result.running:
        if verbose:
            tompytools.generate_message(
                '%s Waiting %i seconds for BLAST' % (refseq_id, wait_seconds))
        time.sleep(wait_seconds)
        if verbose:
            tompytools.generate_message('%s.pollResults()' % refseq_id)
        blast_result.pollResults()
    # check if we found primers
    if verbose:
        tompytools.generate_message('%s.checkSuccess()' % refseq_id)
    blast_result.checkSuccess()
    if verbose:
        tompytools.generate_message('%s.checkSpecificity()' % refseq_id)
    blast_result.checkSpecificity()
    if not (blast_result.noPrimersFound or blast_result.offTargets):
        if verbose:
            tompytools.generate_message('%s.parsePrimers()' % refseq_id)
        blast_result.parsePrimers()
    # return primerBlastResults
    return(blast_result)


# function for looping logic
def run_iterative_primer_blast(
        test_refseq,
        strict_parameters,
        wait_seconds=60,
        verbose=False):
    test_blast_result = run_primer_blast(
        refseq_id=test_refseq,
        blast_parameters=strict_parameters,
        status='strict',
        wait_seconds=wait_seconds,
        verbose=verbose)
    # re-run without intron span
    if test_blast_result.no_intron:
        tompytools.generate_message('Record %s has no introns' % test_refseq)
        test_blast_result.blast_parameters.pop('SPAN_INTRON')
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # re-run with lower clamp requirements
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message('%s: Relaxing GC clamp' % test_refseq)
        test_blast_result.status = 'GC1'
        test_blast_result.blast_parameters['GC_CLAMP'] = '1'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        test_blast_result.status = 'GC0'
        test_blast_result.blast_parameters['GC_CLAMP'] = '0'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # rerun with lower GC requirements
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message('%s: Relaxing GC content' % test_refseq)
        test_blast_result.status = 'GC_content'
        test_blast_result.blast_parameters['PRIMER_MIN_GC'] = '35'
        test_blast_result.blast_parameters['PRIMER_MAX_GC'] = '65'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # rerun with lower primer TM requirements
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message('%s: Relaxing primer TM' % test_refseq)
        test_blast_result.status = 'Low_TM'
        test_blast_result.blast_parameters['PRIMER_MIN_TM'] = '52'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # rerun with lower primer complementarity requirements
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message(
            '%s: Relaxing primer self-complementarity' % test_refseq)
        test_blast_result.status = 'Potential_Dimers'
        test_blast_result.blast_parameters['SELF_ANY'] = '5'
        test_blast_result.blast_parameters['SELF_END'] = '2'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # rerun with near defaults
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message(
            '%s: Using default primer self-complementarity (caution)' %
            test_refseq)
        test_blast_result.status = 'Probable_Dimers'
        test_blast_result.blast_parameters['SELF_ANY'] = '8'
        test_blast_result.blast_parameters['SELF_END'] = '3'
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # rerun without complexity filter
    if test_blast_result.noPrimersFound or test_blast_result.offTargets:
        tompytools.generate_message(
            '%s: Disabling repeat filter' % test_refseq)
        test_blast_result.status = 'No_repeat_filter'
        test_blast_result.blast_parameters.pop('LOW_COMPLEXITY_FILTER')
        test_blast_result = run_primer_blast(
            refseq_id=test_blast_result.RefSeq,
            blast_parameters=test_blast_result.blast_parameters,
            status=test_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    # deal with leftover genes
    if test_blast_result.noPrimersFound:
        test_blast_result.status = 'primer_quality_too_low'
    if test_blast_result.offTargets:
        test_blast_result.status = 'no_specific_primers'
    # return results
    return(test_blast_result)


##############
# PARAMETERS #
##############

# initial blast settings
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

# tell NCBI who we are
my_email = ''
strict_parameters['EMAIL'] = my_email

# try to comply with NCBI usage guidelines:

# Do not poll for any single RID more often than once a minute.
wait_seconds = 60

# Do not contact the server more often than once every three seconds.
max_jobs = int(wait_seconds/3)

# have to increase the recursion limit for joblib
sys.setrecursionlimit(10000)

#############
# TEST CODE #
#############

# gene of interest
# LOC_Os02g41460,XM_015768655,OS02G0623400,G1L3: has introns, off-targets
# LOC_Os05g41760 Os05g0497200, XM_015782903.1, NM_001062476: no introns
test_refseq = 'NM_001062476'
test_iterative_blast_result = run_iterative_primer_blast(
    test_refseq=test_refseq,
    strict_parameters=strict_parameters,
    wait_seconds=10,
    verbose=True)
test_iterative_blast_result_2 = run_iterative_primer_blast(
    test_refseq='XM_015768655',
    strict_parameters=strict_parameters,
    wait_seconds=10,
    verbose=True)



# multiple genes (joblib)
refseq_list = ['NM_001062476', 'XM_015768655']
jobs_to_run = min(len(refseq_list), max_jobs)
blast_results = joblib.Parallel(n_jobs=jobs_to_run, verbose=100)(
    joblib.delayed(run_iterative_primer_blast)(
        test_refseq=x,
        strict_parameters=strict_parameters,
        wait_seconds=wait_seconds,
        verbose=True) for x in refseq_list)

long_gene_list = [
    'AK287830', 'CT837746', 'NM_001066571', 'NM_001061339',
    'NM_001048798', 'NM_001050743', 'CT836522', 'NM_001060130', 'AK243143',
    'NM_001060972', 'NM_001066641', 'NM_001067072', 'NM_001059848']

long_blast_results = joblib.Parallel(n_jobs=10, verbose=100)(
    joblib.delayed(run_iterative_primer_blast)(
        test_refseq=x,
        strict_parameters=strict_parameters,
        wait_seconds=wait_seconds,
        verbose=True) for x in long_gene_list)

# process blast results into CSV etc.
with open('test/test_summary.csv', 'w') as file:
    file.write('RefSeq,status,F,TM_F,R,TM_R,'
               'ProductSize,IntronSize\n')
    for result in long_blast_results:
        file.write(result.csvLine() + '\n')


# some genes eh? NM_001059848