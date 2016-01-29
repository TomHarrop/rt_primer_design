# -*- coding: utf-8 -*-

import requests
from bs4 import BeautifulSoup
import re

primerBlastUrl = 'http://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'

#########
# Class #
#########

# The __init__ and pollResults methods both retrieve the current results page
# from the BLAST server and convert it to html for parsing. The main difference
# is that pollResults prints the html to a file.
#
# Instead of running the methods from the main script and storing them in a
# dict there, they could be called from the main script and stored as
# attributes of the primerBlastResults instance. This would make the
# primerBlastResults instance act like a big dictionary (one:many). To acheive
# this, the __init__ method has to initiate empty paramaters, which are then
# filled out by the subsequent methods (so the methods update the instance
# rather than returning values)


class primerBlastResults:
    '''For retrieving and parsing results from the NCBI Primer-BLAST server.'''
    def __init__(self, job_key, LOC, RefSeq,
                 html=None, url=None, running=None, exceptions=None,
                 noPrimersFound=None, offTargets=None, finalStatus=None,
                 F=None, R=None, TM_F=None, TM_R=None, ProductSize=None,
                 IntronSize=None,
                 blastUrl='http://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'):
        '''Initiliase an instance of primerBlastResults by downloading the html
        for LOC from blastUrl using job_key.
        '''
        self.job_key = job_key
        self.LOC = LOC
        self.blastUrl = blastUrl
        self.RefSeq = RefSeq
        statusPageResponse = requests.get(self.blastUrl,
                                          params={'job_key': self.job_key})
        self.url = statusPageResponse.url
        self.html = BeautifulSoup(statusPageResponse.content, 'lxml')

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

    def printFile(self, subdir):
        '''(primerBlastResults) -> NoneType

        Output html to file named LOC.html

        '''
        with open(subdir + "/" + self.LOC + ".html", 'w') as file:
            print(self.html, file=file)

    def pollResults(self):
        '''(primerBlastResults) -> NoneType

        Retrieve the current status or results page for primerBlastResults and
        update self.html. Output the html to file LOC.html

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
            self.exceptions = 'Exception' in self.html.find(class_='error').text
        else:
            self.exceptions = False

    def checkSpecificity(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.offTargets to true if the retrieved
        primers may not be template-specific.

        '''
        if self.html.find(class_='paramSummary'):
            self.offTargets = 'may not be specific' in self.html.find(class_='paramSummary').text
        else:
            self.offTargets = False

    def checkSuccess(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and update self.noPrimersFound to True if there was a
        warning running primer-BLAST (usually meaning there were no primers
        found for the criteria given, and the criteria should be relaxed)

        '''
        if self.html.find(class_='warning'):
            self.noPrimersFound = 'loosen the selection criteria' in self.html.find(class_='warning').text
        else:
            self.noPrimersFound = False

    def parsePrimers(self):
        '''(primerBlastResults) -> NoneType

        Parse self.html and add results to self.F, self.R, self.TM_F,
        self.TM_R, self.ProductSize and self.IntronSize

        '''
        pF = re.compile('Forward primer\s*\d*\s*([ACTG]*)')
        pR = re.compile('Reverse primer\s*\d*\s*([ACTG]*)')
        pProductSize = re.compile('product length = (\d*)')
        pIntronSize = re.compile('Total intron size(\d*)')
        pTMF = re.compile('Forward primer.*?(\d\d\\.\d\d)')
        pTMR = re.compile('Reverse primer.*?(\d\d\\.\d\d)')
        primerText = self.html.find(class_='prPairDtl').text
        primerTable = self.html.find(class_='prPairInfo').table.text
        self.F = pF.search(primerText).group(1)
        self.R = pR.search(primerText).group(1)
        self.ProductSize = pProductSize.search(primerText).group(1)
        self.IntronSize = pIntronSize.search(primerTable).group(1)
        self.TM_F = pTMF.search(primerTable).group(1)
        self.TM_R = pTMR.search(primerTable).group(1)

    def csvLine(self):
        '''(primerBlastResults) -> str

        Return a line of csv output.        

        '''
        return '{0},{1},{2},{3},{4},{5},{6},{7},{8}'.format(self.LOC, self.RefSeq, self.finalStatus, self.F, self.TM_F, self.R, self.TM_R, self.ProductSize, self.IntronSize)

    def replaceCssLinks(self):
        '''(primerBlastResults) -> NoneType

        Replace the local links to css files in self.html with global ones.

        '''
        cssLinks = self.html.findAll(href=re.compile('css'))
        for link in cssLinks:
            link['href'] = 'http://www.ncbi.nlm.nih.gov/tools/primer-blast/' + link['href']

#############
# Functions #
#############


def runBlast(MsuRefSeq, parameters,
             primerBlastUrl='http://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'):
    '''(dict of LOC:RefSeq, dict of parameters) -> dict of LOC:primerBlastResults

    Submit primer-BLAST jobs to primerBlastUrl for RefSeq IDs in MsuRefSeq
    using parameters. Poll the BLAST server until jobs finish. Update and parse
    the BLAST results. Return dict of LOC:populated primerBlastResults
    instances.

    '''
    # Function for submitting BLAST job. Currently getting a lot of responses
    # that don't contain a job_key. Have to investigate manually.

    def submitBLASTjob(RefSeqID, parameters,
                       primerBlastUrl='http://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'):
        '''(str, dict, str) -> Response

        Post the BLAST job for RefSeqID to primerBlastUrl with parameters and
        return the http Response.

        '''
        post_params = parameters.copy()
        post_params['INPUT_SEQUENCE'] = RefSeqID
        return(requests.get(primerBlastUrl, params=post_params))

    # Function for getting job_key
    def get_job_key(Response):
        '''(str) -> str

        Parse the html in Response and return the job_key

        '''
        soup = BeautifulSoup(Response.content, 'xml')
        if soup.find(NAME='job_key'):
            return(soup.find(NAME='job_key')['VALUE'])

    # Code for runBlast starts here

    from progressbar import ProgressBar, Timer
    import requests
    from bs4 import BeautifulSoup
    import time

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
        # for testing
        # for i in range(10):
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
