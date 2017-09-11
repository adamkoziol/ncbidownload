#!/usr/bin/env python 3
from accessoryFunctions.accessoryFunctions import *
from Bio import Entrez, SeqIO
from time import sleep
__author__ = 'adamkoziol'

Entrez.email = 'adam.koziol@inspection.gc.ca'


class Ncbiquery(object):

    def download(self):
        """
        Download Genbank files corresponding to the supplied accession number from NCBI
        """
        printtime('Downloading and formatting Genbank files from NCBI', self.start)
        from threading import Thread
        # Create and start threads
        for accession in self.accessions:
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.downloadthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
            metadata = MetadataObject()
            metadata.accession = accession
            self.metadata.append(metadata)

        for sample in self.metadata:
            # Create the file name
            sample.genbankfile = os.path.join(self.sequencepath, sample.accession + '.gbk')
            try:
                # Add the sample to the queue
                self.queue.put(sample)
            except (KeyboardInterrupt, SystemExit):
                printtime('Keyboard interrupt. Exiting', self.start)
                quit()
        # Join the threads
        self.queue.join()
        # Parse the Genbank files
        self.parse()

    def downloadthreads(self):
        while True:
            sample = self.queue.get()
            # Attempt to download up to ten times
            for i in range(10):
                try:
                    _ = os.stat(sample.genbankfile).st_size
                    zero = False
                except FileNotFoundError:
                    zero = True
                if zero or not os.path.isfile(sample.genbankfile):
                    # https://stackoverflow.com/a/2083996
                    while True:
                        try:
                            # Use the efetch utility to download the Genbank record in text format
                            # from the nucleotide database
                            handle = Entrez.efetch(db="assembly",
                                                   id=sample.accession,
                                                   rettype="gb",
                                                   retmode="text")
                            # Write the record to file to keep from having to re-download the file if the
                            # script needs to be run multiple times
                            with open(sample.genbankfile, 'w') as genbankfile:
                                genbankfile.write(handle.read())
                            # Sleep in order to keep from overloading NCBI servers with too many requests
                            sleep(0.5)
                        except Exception:
                            continue
                        break
            self.queue.task_done()

    def download_biosample(self):
        """
        Download Genbank files corresponding to the supplied accession number from NCBI
        """
        printtime('Downloading and formatting Genbank files from NCBI', self.start)
        # Make the MetadataObjects with the accession as the sample.accession
        for accession in self.accessions:
            metadata = MetadataObject()
            metadata.accession = accession
            self.metadata.append(metadata)
        for sample in self.metadata:
            # Create the file name
            sample.genbankfile = os.path.join(self.sequencepath, sample.accession + '.gbk')
            if not os.path.isfile(sample.genbankfile):
                # Try to download up to ten times
                for i in range(10):
                    # Determine if the file is empty
                    try:
                        _ = os.stat(sample.genbankfile).st_size
                        zero = False
                    except FileNotFoundError:
                        zero = True
                    if zero or not os.path.isfile(sample.genbankfile):
                        # Method to retry the download up to ten times - call download() once. If an exception occurs,
                        # run the function from https://stackoverflow.com/a/39249214
                        def download():
                            try:
                                # https://www.biostars.org/p/141581/
                                record = Entrez.read(Entrez.esearch(db='assembly',
                                                                    term=sample.accession + '.1[ASAC]',
                                                                    retmode='text'))
                                handle = str()
                                for id_list in record['IdList']:
                                    # Get Assembly Summary
                                    esummary_record = Entrez.read(Entrez.esummary(db="assembly",
                                                                                  id=id_list,
                                                                                  report="full"))
                                    # Parse Accession Name from Summary
                                    biosample = \
                                        esummary_record['DocumentSummarySet']['DocumentSummary'][0]['BioSampleId']
                                    # Fetch record from BioSample db
                                    handle = Entrez.efetch(db="biosample",
                                                           id=biosample,
                                                           rettype='gb')
                                # Write the record to file to keep from having to re-download the file if the script
                                # needs to be run multiple times
                                with open(sample.genbankfile, 'w') as genbankfile:
                                    genbankfile.write(handle.read())
                            except Exception:
                                # Sleep to avoid overloading the NCBI servers
                                sleep(1)
                                # Try to download again
                                download()
                        # Run the download function
                        download()
                    else:
                        break
        # Parse the BioSample record
        self.parsebiosample()

    def parse(self):
        """
        Parse the Genbank files, and extract the 'host' feature from strains containing that feature
        """
        printtime('Parsing Genbank files', self.start)
        # Open the file of accessions and corresponding hosts
        with open(os.path.join(self.path, 'country_full.txt'), 'w') as hosts:
            for sample in self.metadata:
                # Read in the Genbank record from disk
                record = SeqIO.read(sample.genbankfile, 'genbank')
                # Iterate through the features in the Genbank file
                for feature in record.features:
                    # Find the 'host' key in the features.qualifiers dictionary
                    try:
                        # host = feature.qualifiers['host'][0]
                        # source = feature.qualifiers['isolation_source'][0]
                        country = feature.qualifiers['country'][0]
                        # Write the hosts to file
                        hosts.write('{}\t{}\n'.format(sample.accession, country))
                    # Sometimes the host is not recorded. Pass
                    except KeyError:
                        pass

    def parsebiosample(self):
        """
        Parse the Genbank files, and extract the 'host' feature from strains containing that feature
        """
        import xml.etree.ElementTree as ElementTree
        printtime('Parsing Genbank files', self.start)
        # Open the file of accessions and corresponding hosts
        with open(os.path.join(self.path, 'accession_organism_host.txt'), 'w') as hosts:
            for sample in self.metadata:
                # Initialise attributes with default values to be overwritten if
                sample.host = 'missing'
                sample.organism = 'unknown'
                # Load the XML file using ElementTree
                xml = ElementTree.parse(sample.genbankfile)
                # Get the root of the XML tree for parsing the file
                root_element = xml.getroot()
                # Iterate through the tree - the first child for these files will be 'BioSampleSet'
                for biosample in root_element:
                    # Find the headers within 'BioSampleSet
                    for header in biosample:
                        # There are two headers containing data of interest: Description and Attributes
                        if header.tag == 'Description':
                            # Iterate through the description sub-headers
                            for description in header:
                                # Find the description tag with the 'Organism' keyword
                                if description.tag == 'Organism':
                                    # Set the organism name as the 'taxonomy_name' attribute
                                    sample.organism = description.attrib['taxonomy_name']
                        elif header.tag == 'Attributes':
                            # Iterate through all the attribute sub-headers
                            for attribute in header:
                                # Find the attribute matching 'specific host'
                                if attribute.attrib['attribute_name'] == 'specific host':
                                    # Set the host as the text field from the 'specific host' attribute
                                    sample.host = attribute.text
                                # Hosts are sometimes specified only by their taxonomy ID
                                elif attribute.attrib['attribute_name'] == 'Host taxon ID':
                                    # Extract the id
                                    taxid = attribute.text
                                    # Fetch the taxonomy information from NCBI
                                    taxonomy = Entrez.read(Entrez.efetch(db='Taxonomy',
                                                                         id=taxid))
                                    # Set the host
                                    sample.host = taxonomy[0]['ScientificName']
                                # Another method hosts can be specified
                                elif attribute.attrib['attribute_name'] == 'host':
                                    sample.host = attribute.text
                hosts.write('{}\t{}\t{}\n'.format(sample.accession, sample.organism, sample.host))

    def __init__(self, args):
        from queue import Queue
        self.path = args.path
        self.sequencepath = os.path.join(args.sequencepath, '')
        make_path(self.sequencepath)
        self.file = os.path.join(self.path, args.file)
        self.assembly = args.assembly
        # Read in all the accessions from the file
        try:
            with open(self.file, 'r') as accessions:
                self.accessions = accessions.read().splitlines()
        except IOError:
            print('Cannot find file of accessions as provided: {}'.format(self.file))
            quit()
        self.start = args.start
        self.queue = Queue(maxsize=5)
        self.metadata = list()
        if not self.assembly:
            self.download()
        else:
            self.download_biosample()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    # Parser for arguments
    parser = ArgumentParser(description='Download Genbank files from NCBI to search for the host of bacterial strains')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to store downloaded Genbank files')
    parser.add_argument('-f', '--file',
                        help='Name of file of accession numbers')
    parser.add_argument('-a', '--assembly',
                        action='store_true',
                        help='If the input accessions are for the "assembly" database, different functions must be used'
                        )
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()
    # Run it
    Ncbiquery(arguments)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')
