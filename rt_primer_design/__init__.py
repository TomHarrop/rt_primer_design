#!/usr/bin/env python3
# -*- coding: utf-8 -*-


################
# DEPENDENCIES #
################

from bs4 import BeautifulSoup
import joblib
import re
import requests
import sys
import time
import tompytools


############
# PROVIDES #
############

class PrimerBlastResult:
    """Submit BLAST query to NCBI and parse results.

    Attributes:
        F: Forward primer
        R: Reverse primer
        TM_F: TM of the forward primer
        TM_R: TM of the reverse primer
        exceptions: True if the string 'Exception' is in the text of 'error'
            class in html
        html: result of GET parsed with BeautifulSoup lxml parser
        intron_size: expected size of intron(s) in the PCR product from genomic
            DNA
        job_key: job_key for retrieving results from NCBI server.
        no_intron: True if the string 'junction cannot be found' is in the text
            of 'info' class in the html
        no_primers_found: True if 'No primers were found' is in the text of the
            'info' class in html or if 'loosen the selection criteria' is in
            the 'warning' class of the html
        product_size: expected PCR product size
        off_targets: True if string 'may not be specific' is in the text of the
            'paramSummary' class in the html
        running: True if the job is still running on the NCBI server
        url: complete URL used for the GET query by submit_blast_request
        user_seqloc: USER_SEQLOC values for similar templates according to the
            NCBI server

    """
    def __init__(
        self, ref_seq, status, blast_parameters,
        F=None,
        R=None,
        TM_F=None,
        TM_R=None,
        exceptions=None,
        html=None,
        intron_size=None,
        job_key=None,
        no_intron=None,
        no_primers_found=None,
        off_targets=None,
        product_size=None,
        running=None,
        url=None,
        user_seqloc=None,
        blast_url=('https://www.ncbi.nlm.nih.gov/tools/'
                   'primer-blast/primertool.cgi')):
        """Init PrimerBlastResult.

        Submits the BLAST query for ref_seq to blast_url using blast_parameters
        and returns a PrimerBlastResult with attributes blast_url, ref_seq,
        blast_parameters, status, url, html, job_key, running.

        Args:
            ref_seq: RefSeq sequence ID to design primers for.
            status: the status of the primer set, e.g. 'strict' if submitting
                the initial BLAST search with strict parameters.
            blast_parameters: the parameters to use for designing the primers.
        """
        # initiate the object
        self.blast_url = blast_url
        self.ref_seq = ref_seq
        self.status = status
        self.blast_parameters = blast_parameters.copy()
        self.blast_parameters['INPUT_SEQUENCE'] = self.ref_seq

        # submit the BLAST request
        self.submit_blast_request()

    def __eq__(self, other):
        """Compare PrimerBlastResult objects by job_key.

        Args:
            other: the PrimerBlastResult to compare to self.

        Returns:
            A boolean value indicating if the job_keys are identical.
        """
        return self.job_key == other.job_key

    def __str__(self):
        """Print formatted html result."""
        print(self.html)

    def submit_blast_request(self):
        """Run the BLAST query on NCBI primer-blast.

        Submits the BLAST query then calls get_job_key and poll_results
        methods.
        """
        # submit the BLAST request and get the response page
        blast_result = requests.get(
            self.blast_url,
            params=self.blast_parameters)
        self.url = blast_result.url
        self.html = BeautifulSoup(blast_result.content, 'lxml')

        # get a job_key
        self.get_job_key()

        # poll the results once to get the finished page
        self.poll_results()

    def get_job_key(self):
        """Retrieve the job_key and store in self.job_key"""
        # first choice is the proper "job_key" tag
        if self.html(attrs={'name': 'job_key'}):
            self.job_key = self.html.find(attrs={'name': 'job_key'})['value']
            # print('%s: Found job_key %s in attrs' %
            #       (self.ref_seq, self.job_key))

        # otherwise, try to parse the 'Job id' from `breadcrumb` with regex :(
        if self.html.find(id='breadcrumb'):
            bc_text = self.html.find(id='breadcrumb').text
            bc_search = re.compile(r'Job id\=(\S+).*')
            if bc_search.search(bc_text):
                self.job_key = bc_search.search(bc_text).group(1)
                # print('%s: Parsed job_key %s from breadcrumb' %
                #       (self.ref_seq, self.job_key))

        # otherwise, give up and submit the blast search again
        if not hasattr(self, "job_key"):
            print('''%s: Couldn't parse job_key''' % self.ref_seq)
            print(self.html)
            print('%s: Waiting 60 seconds and '
                  're-submitting the BLAST request' % self.ref_seq)
            time.sleep(60)
            self.submit_blast_request()

    def poll_results(self):
        """Retrieve the current html from NCBI and update self.html and
        self.running.
        """
        status_page_response = requests.get(
            self.blast_url,
            params={'job_key': self.job_key})
        self.html = BeautifulSoup(status_page_response.content, 'lxml')
        self.check_running()

    def check_running(self):
        """Parse self.html and update self.running."""
        if self.html.find(class_='odd'):
            self.running = 'Running' in self.html.find(class_='odd').text
        else:
            self.running = False

    def check_introns(self):
        """Parse self.html and update self.no_intron."""
        if self.html.find(class_='info'):
            self.no_intron = ('junction cannot be found' in
                              self.html.find(class_='info').text)
        else:
            self.no_intron = False

    def check_specificity(self):
        """Parse self.html and update self.off_targets."""
        if self.html.find(class_='paramSummary'):
            self.off_targets = ('may not be specific' in
                                self.html.find(class_='paramSummary').text)
        else:
            self.off_targets = False

    def check_success(self):
        """Parse self.html and update self.no_primers_found.

        Parse self.html and update self.no_primers_found based on the 'warning'
        and 'info' html classes. This usually means no primers were found for
        the given blast_parameters, which should be relaxed.
        """
        warning_list = []
        if self.html.find(class_='warning'):
            warning_list.append(
                ('loosen the selection criteria' in
                 self.html.find(class_='warning').text))

        if self.html.find(class_='info'):
            warning_list.append(
                ('No primers were found' in
                 self.html.find(class_='info').text))

        if any(x for x in warning_list):
            self.no_primers_found = True
        else:
            self.no_primers_found = False

    def parse_primers(self):
        """Parse self.html for primer details."""
        primer_table = self.html.find(class_='prPairInfo').table.text

        F_search = re.compile('Forward primer\s*\d*\s*([ACTG]+)')
        R_search = re.compile('Reverse primer\s*\d*\s*([ACTG]+)')
        product_size_search = re.compile('Product length(\d*)')
        intron_size_search = re.compile('Total intron size(\d*)')
        TM_F_search = re.compile('Forward primer.*?(\d\d\\.\d\d)')
        TM_R_search = re.compile('Reverse primer.*?(\d\d\\.\d\d)')

        self.F = F_search.search(primer_table).group(1)
        self.R = R_search.search(primer_table).group(1)
        self.product_size = product_size_search.search(primer_table).group(1)
        self.intron_size = intron_size_search.search(primer_table).group(1)
        self.TM_F = TM_F_search.search(primer_table).group(1)
        self.TM_R = TM_R_search.search(primer_table).group(1)

    def csv_line(self):
        """Print a line of csv.

        Use the primer attributes found by parse_primers to print a line of CSV
        for writing to file.

        Returns:
            A str of text in csv format.
        """
        if not hasattr(self, 'F'):
            return '{0},{1},,,,,,'.format(self.ref_seq, self.status)

        return '{0},{1},{2},{3},{4},{5},{6},{7}'.format(
            self.ref_seq,       # 0
            self.status,        # 1
            self.F,             # 2
            self.TM_F,          # 3
            self.R,             # 4
            self.TM_R,          # 5
            self.product_size,  # 6
            self.intron_size)   # 7

    def replace_css_links(self):
        """Replace the local css links in self.html with absolute paths."""
        css_links = self.html.findAll(href=re.compile('css'))
        for link in css_links:
            link['href'] = ('https://www.ncbi.nlm.nih.gov/'
                            'tools/primer-blast/' + link['href'])

    def check_similar_templates(self):
        """Check for similar templates and re-submit BLAST if necessary.

        If NCBI reports that the template is "highly similar" to another
        sequence, parse the RefSeq ID of that sequence and resubmit the BLAST
        query with the extra parameters 'TRY_USER_GUIDE' and 'USER_SEQLOC',
        and update the user_seqloc attribute.
        """
        if self.html.find(id='expl'):
                if (('Your PCR template is highly similar to the following'
                     ' sequence') in self.html.find(id='expl').text):
                    self.user_seqloc = [
                        x['value'] for x in self.html.find_all(
                            name='input',
                            type='checkbox',
                            attrs={'name': 'USER_SEQLOC'})]
                    print("%s: Found similar sequences: %s" %
                          (self.ref_seq, self.user_seqloc))
                    self.blast_parameters['TRY_USER_GUIDE'] = 'yes'
                    self.blast_parameters['USER_SEQLOC'] = self.user_seqloc

                    # this is a new BLAST request
                    self.submit_blast_request()

    def print_file(self, output_subdirectory):
        """Print html to a file in output_subdirectory"""
        output_file = output_subdirectory + "/" + self.ref_seq + ".html"
        with open(output_file, 'w') as file:
            print(self.html, file=file)

def run_primer_blast(
        ref_seq,
        blast_parameters,
        status,
        wait_seconds=60,
        verbose=False):
    """Run NCBI primer-blast and wait for results.

    Submits a BLAST query for ref_seq using blast_parameters and an initial
    status. Waits for BLAST server to finish query, polling every wait_seconds.
    wait_seconds should be 60 to as per NCBI usage guidelines. Checks results
    for exon/exon junction. Checks for similar templates and waits for new
    BLAST query to finish if necessary. Lastly, checks whether primers were
    found and if they are specific. If so, parses the primers.

    Args:
        ref_seq: RefSeq sequence ID to design primers for.
        status: the status of the primer set, e.g. 'strict' if submitting
            the initial BLAST search with strict parameters.
        blast_parameters: the parameters to use for designing the primers.
        wait_seconds: how long to wait between polling the BLAST server.
            NCBI usage guidelines state "Do not poll for any single RID more
            often than once a minute".
        verbose: If True, print stepwise status messages.

    Returns:
        A PrimerBlastResult with attributes no_intron, no_primers_found and
        off_targets.
    """
    if verbose:
        tompytools.generate_message('Running BLAST for %s' % ref_seq)
        print(
            'ref_seq=%s\t\t\t'
            'blast_parameters=%s\t\t\t'
            'status=%s' %
            (ref_seq, blast_parameters, status))

    # run the BLAST search
    blast_result = PrimerBlastResult(
        ref_seq=ref_seq,
        blast_parameters=blast_parameters,
        status=status)

    # wait for job to finish
    while blast_result.running:
        if verbose:
            tompytools.generate_message(
                '%s: Waiting %i seconds for BLAST' % (ref_seq, wait_seconds))
        time.sleep(wait_seconds)
        if verbose:
            tompytools.generate_message('%s.poll_results()' % ref_seq)
        blast_result.poll_results()

    # check for exon/exon junction
    if verbose:
        tompytools.generate_message('%s.check_introns()' % ref_seq)
    blast_result.check_introns()

    # check for similar sequences and wait for re-run to finish if required
    if verbose:
        tompytools.generate_message('%s.check_similar_templates()' % ref_seq)
    blast_result.check_similar_templates()
    while blast_result.running:
        if verbose:
            tompytools.generate_message(
                '%s: Waiting %i seconds for BLAST' % (ref_seq, wait_seconds))
        time.sleep(wait_seconds)
        if verbose:
            tompytools.generate_message('%s.poll_results()' % ref_seq)
        blast_result.poll_results()

    # check if we found primers
    if verbose:
        tompytools.generate_message('%s.check_success()' % ref_seq)
    blast_result.check_success()
    if verbose:
        tompytools.generate_message('%s.check_specificity()' % ref_seq)
    blast_result.check_specificity()
    if not (blast_result.no_primers_found or blast_result.off_targets):
        if verbose:
            tompytools.generate_message('%s.parse_primers()' % ref_seq)
        blast_result.parse_primers()

    # finished
    return blast_result


def iterate_primer_blast(
        ref_seq,
        starting_parameters,
        wait_seconds=60,
        verbose=False):
    """Progressively relax BLAST parameters until primers are found.

    Runs a BLAST query for ref_seq with starting_parameters using
    run_primer_blast(). Checks resuts and re-runs with relaxed parameters if
    the gene has no introns, no primers were found, or the primers are not
    specific.

    Args:
        ref_seq: RefSeq sequence ID to design primers for.
        starting_parameters: the strict parameters, which will be progressively
            relaxed until primers are found.
        wait_seconds: how long to wait between polling the BLAST server.
            NCBI usage guidelines state "Do not poll for any single RID more
            often than once a minute".
        verbose: If True, print stepwise status messages.
    Returns:
        A PrimerBlastResult with status and parsed primers.
    """
    # start iteration with starting_parameters
    iterative_blast_result = run_primer_blast(
        ref_seq=ref_seq,
        blast_parameters=starting_parameters,
        status='strict',
        wait_seconds=wait_seconds,
        verbose=verbose)

    # check for intron, pop and re-run if necessary
    if iterative_blast_result.no_intron:
        tompytools.generate_message('Record %s has no introns' % ref_seq)
        iterative_blast_result.blast_parameters.pop('SPAN_INTRON')
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # re-run with lower clamp requirements
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message('%s: Relaxing GC clamp' % ref_seq)
        iterative_blast_result.status = 'GC1'
        iterative_blast_result.blast_parameters['GC_CLAMP'] = '1'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        iterative_blast_result.status = 'GC0'
        iterative_blast_result.blast_parameters['GC_CLAMP'] = '0'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # re-run with lower GC requirements
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message('%s: Relaxing GC content' % ref_seq)
        iterative_blast_result.status = 'GC_content'
        iterative_blast_result.blast_parameters['PRIMER_MIN_GC'] = '35'
        iterative_blast_result.blast_parameters['PRIMER_MAX_GC'] = '65'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # re-run with lower primer TM requirements
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message('%s: Relaxing primer TM' % ref_seq)
        iterative_blast_result.status = 'Low_TM'
        iterative_blast_result.blast_parameters['PRIMER_MIN_GC'] = '35'
        iterative_blast_result.blast_parameters['PRIMER_MAX_GC'] = '65'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # rerun with lower primer complementarity requirements
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message(
            '%s: Relaxing primer self-complementarity' % ref_seq)
        iterative_blast_result.status = 'Potential_Dimers'
        iterative_blast_result.blast_parameters['SELF_ANY'] = '5'
        iterative_blast_result.blast_parameters['SELF_END'] = '2'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # rerun with near defaults
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message(
            '%s: Using default primer self-complementarity (caution)' %
            ref_seq)
        iterative_blast_result.status = 'Probable_Dimers'
        iterative_blast_result.blast_parameters['SELF_ANY'] = '8'
        iterative_blast_result.blast_parameters['SELF_END'] = '3'
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # rerun without complexity filter
    if (iterative_blast_result.no_primers_found or
            iterative_blast_result.off_targets):
        tompytools.generate_message(
            '%s: Disabling repeat filter' %
            ref_seq)
        iterative_blast_result.status = 'No_repeat_filter'
        iterative_blast_result.blast_parameters.pop('LOW_COMPLEXITY_FILTER')
        iterative_blast_result = run_primer_blast(
            ref_seq=iterative_blast_result.ref_seq,
            blast_parameters=iterative_blast_result.blast_parameters,
            status=iterative_blast_result.status,
            wait_seconds=wait_seconds,
            verbose=verbose)

    # deal with leftover genes
    if iterative_blast_result.no_primers_found:
        iterative_blast_result.status = 'primer_quality_too_low'
    if iterative_blast_result.off_targets:
        iterative_blast_result.status = 'no_specific_primers'

    # finished
    return iterative_blast_result


def multiple_primer_blast(
        ref_seq_list,
        starting_parameters,
        wait_seconds=60,
        verbose=False,
        n_jobs=10):
    """Wrapper to joblib.Parallel for running multiple BLAST queries

    Launches processes  and runs iterate_primer_blast for each gene in
    ref_seq_list using starting_parameters. By default runs 10 jobs
    simultaneously.

    Args:
        ref_seq_list: list of RefSeq IDs
        starting_parameters: the strict parameters, which will be progressively
            relaxed until primers are found.
        wait_seconds: how long to wait between polling the BLAST server.
            NCBI usage guidelines state "Do not poll for any single RID more
            often than once a minute".
        verbose: If True, print stepwise status messages.
        n_jobs: number of simultaneous jobs. NCBI usage guidelines state "Do
            not contact the server more often than once every three seconds",
            so n_jobs should be <= 20 for and wait_seconds should be >= 60
            for long jobs.

    Returns:
        A list of PrimerBlastResult objects, one for each gene in
        ref_seq_list.
    """
    if verbose:
        print("Using system recursion limit: %i" % sys.getrecursionlimit())

    joblib_verbose = 100 if verbose else 0
    jobs_to_run = min(n_jobs, len(ref_seq_list))

    parallel_blast_results = joblib.Parallel(
        n_jobs=jobs_to_run,
        verbose=joblib_verbose)(
        joblib.delayed(iterate_primer_blast)(
            ref_seq=x,
            starting_parameters=starting_parameters,
            wait_seconds=wait_seconds,
            verbose=verbose) for x in ref_seq_list)
    return parallel_blast_results
