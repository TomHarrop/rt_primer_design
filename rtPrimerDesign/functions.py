#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import requests
import re
import time
from bs4 import BeautifulSoup
from progressbar import ProgressBar, Timer


#############
# Variables #
#############

primerBlastUrl = ('https://www.ncbi.nlm.nih.gov/tools/'
                  'primer-blast/primertool.cgi')


#########
# Class #
#########

# The __init__ and pollResults methods both retrieve the current results page
# from the BLAST server and convert it to html for parsing.

# The methods are called from the main script and the results are stored as
# attributes of the primerBlastResults instance. This makes the
# primerBlastResults instance act like a big dictionary (one:many). To acheive
# this, the __init__ method has to initiate empty paramaters, which are then
# filled out by the methods (so the methods update the instance rather than
# returning values)

class primerBlastResults:
    '''For retrieving and parsing results from the NCBI Primer-BLAST server.'''
    def __init__(self, RefSeq, blast_parameters,
                 job_key=None,
                 html=None,
                 url=None,
                 running=None,
                 exceptions=None,
                 noPrimersFound=None,
                 offTargets=None,
                 finalStatus=None,
                 F=None,
                 R=None,
                 TM_F=None,
                 TM_R=None,
                 ProductSize=None,
                 IntronSize=None,
                 user_seqloc=None,
                 blastUrl=primerBlastUrl):

        '''
        Initiliase an instance of primerBlastResults by submitting the BLAST
        search and downloading the html response
        '''
        self.blastUrl = blastUrl
        self.RefSeq = RefSeq
        self.blast_parameters = blast_parameters.copy()

        # submit the BLAST request and get the response page
        self.blast_parameters['INPUT_SEQUENCE'] = self.RefSeq
        blast_result = requests.get(
            self.blastUrl,
            params=self.blast_parameters)
        self.url = blast_result.url
        self.html = BeautifulSoup(blast_result.content, 'lxml')

        # get a job_key
        self.get_job_key()

        # poll the results once to get the finished page
        self.pollResults()

    def __eq__(self, other):
        '''(primerBlastResults, primerBlastResults) -> bool

        Return whether primerBlastResult has the same job_key as other

        '''
        return self.job_key == other.job_key

    def __str__(self):
        '''(primerBlastResults) -> str

        Return the html stored in primerBlastResults.html

        '''
        return self.html

    def get_job_key(self):
        '''(primerBlastResults) -> NoneType

        Retrieve the job_key and store in self.job_key

        '''

        # get the job_key. first choice is the proper "job_key" tag
        if self.html(attrs={'name': 'job_key'}):
            self.job_key = self.html.find(attrs={'name': 'job_key'})['value']
        # otherwise, try to parse the 'Job id' from `breadcrumb` with regex :(
        elif self.html.find(id='breadcrumb'):
            bc_strings = [
                x for x in self.html.find(id='breadcrumb').stripped_strings]
            bc_search = re.compile(r'Job id\=(\S+)$')
            for bc_string in bc_strings:
                if bc_search.search(bc_string):
                    self.job_key = bc_search.search(bc_string).groups()[0]

    def printFile(self, subdir):
        '''(primerBlastResults) -> NoneType

        Output html to file named LOC.html

        '''
        with open(subdir + "/" + self.LOC + ".html", 'w') as file:
            print(self.html, file=file)

    def pollResults(self):
        '''(primerBlastResults) -> NoneType

        Retrieve the current status or results page for primerBlastResults and
        update self.html.

        '''
        statusPageResponse = requests.get(self.blastUrl,
                                          params={'job_key': self.job_key})
        self.html = BeautifulSoup(statusPageResponse.content, 'lxml')

    def checkRunning(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.running to True if the job is still
        running

        '''
        if self.html.find(class_='odd'):
            self.running = 'Running' in self.html.find(class_='odd').text
        else:
            self.running = False

    def checkExceptions(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.exceptions to True if there was an
        error running primer-BLAST (usually because the gene has no introns)

        '''
        if self.html.find(class_='error'):
            self.exceptions = ('Exception' in
                               self.html.find(class_='error').text)
        elif self.html.find(class_='info'):
            self.exceptions = ('junction cannot be found' in
                               self.html.find(class_='info').text)
        else:
            self.exceptions = False

    def checkSpecificity(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.offTargets to true if the retrieved
        primers may not be template-specific.

        '''
        if self.html.find(class_='paramSummary'):
            self.offTargets = ('may not be specific' in
                               self.html.find(class_='paramSummary').text)
        else:
            self.offTargets = False

    def checkSuccess(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.noPrimersFound to True if there was a
        warning running primer-BLAST (usually meaning there were no primers
        found for the criteria given, and the criteria should be relaxed)

        '''
        if self.html.find(class_='warning'):
            self.noPrimersFound = ('loosen the selection criteria' in
                                   self.html.find(class_='warning').text)
        elif self.html.find(class_='info'):
            self.noPrimersFound = ('No primers were found' in
                                   self.html.find(class_='info').text)
        else:
            self.noPrimersFound = False

    def parsePrimers(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and add results to self.F, self.R, self.TM_F,
        self.TM_R, self.ProductSize and self.IntronSize

        '''
        pF = re.compile('Forward primer\s*\d*\s*([ACTG]+)')
        pR = re.compile('Reverse primer\s*\d*\s*([ACTG]+)')
        pProductSize = re.compile('Product length(\d*)')
        pIntronSize = re.compile('Total intron size(\d*)')
        pTMF = re.compile('Forward primer.*?(\d\d\\.\d\d)')
        pTMR = re.compile('Reverse primer.*?(\d\d\\.\d\d)')
        primerTable = self.html.find(class_='prPairInfo').table.text
        self.F = pF.search(primerTable).group(1)
        self.R = pR.search(primerTable).group(1)
        self.ProductSize = pProductSize.search(primerTable).group(1)
        self.IntronSize = pIntronSize.search(primerTable).group(1)
        self.TM_F = pTMF.search(primerTable).group(1)
        self.TM_R = pTMR.search(primerTable).group(1)

    def csvLine(self):
        '''(primerBlastResults) -> str

        Return a line of csv output.

        '''
        return '{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(
            self.LOC, self.RefSeq, self.finalStatus, self.F, self.TM_F, self.R,
            self.TM_R, self.ProductSize, self.IntronSize)

    def replaceCssLinks(self):
        '''(primerBlastResults) -> NoneType

        Replace the local links to css files in self.html with global ones.

        '''
        cssLinks = self.html.findAll(href=re.compile('css'))
        for link in cssLinks:
            link['href'] = ('http://www.ncbi.nlm.nih.gov/tools/primer-blast/' +
                            link['href'])

    def check_similar_templates(self):
        '''(primerBlastResults) -> NoneType

        Check for similar templates and store in self.user_seqloc. Rerun
        BLAST with similar sequences enabled.

        '''
        if (('Your PCR template is highly similar '
             'to the following sequence') in self.html.find(id='expl').text):
            self.user_seqloc = [
                x['value'] for x in self.html.find_all(
                    name='input',
                    type='checkbox',
                    attrs={'name': 'USER_SEQLOC'})]
            post_parameters = self.blast_parameters.copy()
            post_parameters['TRY_USER_GUIDE'] = 'yes'
            post_parameters['USER_SEQLOC'] = self.user_seqloc

            # this is effectively a new BLAST search so expect a new job_key
            blast_result = requests.get(
                self.blastUrl,
                params=post_parameters)
            self.html = BeautifulSoup(blast_result.content, 'lxml')
            self.url = blast_result.url
            self.get_job_key()
            self.pollResults()


#############
# Functions #
#############

def run_primer_blast(refseq_ids):
    '''(list) -> dict of refseq_ids:primerBlastResults

    Initialise a primerBlastResults for each entry in refseq_ids. Run the
    BLAST search, poll results until they finish and check for problems.
    The problems will be checked here but deciding what to do with the primers
    should happen in __main__

    PROBLEM                         ACTION
    no exon/exon junction           store exception and process separately
    similar sequences               parse USER_SEQLOC and resubmit query
    no primers found                store in self.noPrimersFound
    primers not specific            store in self.offTargets

    '''
    pass


def submitBLASTjob(RefSeqID, parameters,
                   primerBlastUrl=primerBlastUrl):
    '''(str, dict, str) -> Response

    Post the BLAST job for RefSeqID to primerBlastUrl with parameters and
    return the http Response.

    '''
    post_params = parameters.copy()
    post_params['INPUT_SEQUENCE'] = RefSeqID
    return(requests.get(primerBlastUrl, params=post_params))


def get_job_key(Response):
    '''(str) -> str

    Parse the html in Response and return the job_key

    '''
    soup = BeautifulSoup(Response.content, 'lxml')
    # first choice is the proper "job_key" tag
    if soup.find(attrs={'name': 'job_key'}):
        return(soup.find(attrs={'name': 'job_key'})['value'])
    # otherwise, try to parse the 'Job id' from `breadcrumb` with regex :(
    elif soup.find(id='breadcrumb'):
        bc_strings = [
            x for x in soup.find(id='breadcrumb').stripped_strings]
        bc_search = re.compile(r'Job id\=(\S+)$')
        for bc_string in bc_strings:
            if bc_search.search(bc_string):
                return(bc_search.search(bc_string).groups()[0])


def runBlast(MsuRefSeq, parameters,
             primerBlastUrl=primerBlastUrl):

    '''(dict of LOC:RefSeq, dict of parameters) ->
    dict of LOC:primerBlastResults

    Submit primer-BLAST jobs to primerBlastUrl for RefSeq IDs in MsuRefSeq
    using parameters. Poll the BLAST server until jobs finish. Update and parse
    the BLAST results. Return dict of LOC:populated primerBlastResults
    instances.

    '''

    # Submit BLAST jobs
    print('Submitting BLAST jobs.')
    BlastResponses = {}
    pbar = ProgressBar(maxval=len(MsuRefSeq)).start()
    i = 0
    for key in MsuRefSeq:
        BlastResponses[key] = submitBLASTjob(MsuRefSeq[key], parameters,
                                             primerBlastUrl)
        i = i + 1
        pbar.update(i)
    pbar.finish()

    # Get job_key for each BLAST job
    jobKeys = {}
    for key in BlastResponses:
        jobKeys[key] = get_job_key(BlastResponses[key])
        # print(key)
        # print(jobKeys[key])

    # parse results using class primerBlastResults

    # download the responses
    ResultsPages = {}
    for key in jobKeys:
        ResultsPages[key] = primerBlastResults(
            job_key=jobKeys[key], LOC=key, RefSeq=MsuRefSeq[key])

    # Check if jobs are finished
    for key in ResultsPages:
        ResultsPages[key].checkRunning()

    # Re-poll unfinished jobs until they all finish
    statuses = 0
    for key in ResultsPages:
        statuses += ResultsPages[key].running

    while statuses >= 1:
        print('Waiting 60 seconds for jobs to complete')
        timer = ProgressBar(widgets=[Timer()], maxval=60).start()
        for i in range(60):
            time.sleep(1)
            timer.update(i)
        timer.finish()
        for key in ResultsPages:
            if ResultsPages[key].running:
                ResultsPages[key].pollResults()
                ResultsPages[key].checkRunning()
        statuses = 0
        for key in ResultsPages:
            statuses += ResultsPages[key].running

    # Run the other methods to parse results
    for key in ResultsPages:
        ResultsPages[key].checkExceptions()
        ResultsPages[key].checkSuccess()
        ResultsPages[key].checkSpecificity()

    return ResultsPages
