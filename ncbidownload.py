#!/usr/bin/env python 3
import pandas
from Bio import Entrez
from accessoryfunctions.accessoryFunctions import *
__author__ = 'adamkoziol'

Entrez.email = 'adam.koziol@inspection.gc.ca'


class Download(object):

    def excelparse(self):
        """
        Parses input excel file, and creates objects with headers as keys, and cell data as values for each row
        """
        printtime('Loading accessions from file', self.start)
        # A dictionary to store the parsed excel file in a more readable format
        nesteddictionary = dict()
        # Use pandas to read in the excel file, and subsequently convert the pandas data frame to a dictionary
        # (.to_dict()). Only read the first fourteen columns (parse_cols=range(14)), as later columns are not
        # relevant to this script
        dictionary = pandas.read_excel(self.file).to_dict()
        # Iterate through the dictionary - each header from the excel file
        for header in dictionary:
            # Sample is the primary key, and value is the value of the cell for that primary key + header combination
            for sample, value in dictionary[header].items():
                # Update the dictionary with the new data
                try:
                    nesteddictionary[sample].update({header: value})
                # Create the nested dictionary if it hasn't been created yet
                except KeyError:
                    nesteddictionary[sample] = dict()
                    nesteddictionary[sample].update({header: value})
        # Create objects for each of the samples, rather than using a nested dictionary. It may have been possible to
        # skip the creation of the nested dictionary, and create the objects from the original dictionary, but there
        # seemed to be too many possible places for something to go wrong
        for line in nesteddictionary:
            # Create an object for each sample
            metadata = MetadataObject()
            # Set the name of the metadata to be the primary key for the sample from the excel file
            metadata.name = line
            # Find the headers and values for every sample
            for header, value in nesteddictionary[line].items():
                # Create each attribute - use the header (in lowercase, and spaces removed) as the attribute name,
                # and the value as the attribute value
                setattr(metadata, header.replace(' ', '').lower(), value)
            # Append the object to the list of objects
            self.metadata.append(metadata)
        self.download()

    def download(self):
        """
        Download Genbank files corresponding to the supplied accession number from NCBI
        """
        printtime('Downloading and formatting Genbank files from NCBI', self.start)
        from threading import Thread
        # Create and start threads
        for _ in self.metadata:
            # Send the threads to the appropriate destination function
            threads = Thread(target=self.downloadthreads, args=())
            # Set the daemon to true - something to do with thread management
            threads.setDaemon(True)
            # Start the threading
            threads.start()
        for sample in self.metadata:
            # Create the file name
            sample.genbankfile = os.path.join(self.genbankpath, sample.accession + '.gbk')
            # Add the sample to the queue
            self.queue.put(sample)
        # Join the threads
        self.queue.join()
        # Parse the Genbank files
        self.fastaparse()

    def downloadthreads(self):
        from time import sleep
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

    def fastaparse(self):
        """
        Parse the Genbank files to extract the desired FASTA sequences in the correct orientation
        """
        from Bio import SeqIO
        printtime('Parsing Genbank files', self.start)
        for sample in self.metadata:
            # Create the name of the FASTA-formatted output file
            sample.outputfile = os.path.join(self.sequencepath, '{}_{}.fasta'.format(sample.gene, sample.accession))
            # Read in the Genbank record from disk
            record = SeqIO.read(sample.genbankfile, 'genbank')
            # Set the header to be simply 'gene_accession'
            record.id = '{}_{}'.format(sample.gene, sample.accession)
            record.name = ''
            record.description = ''
            # Extract the sequence desired from the whole sequence (as the start position provided is 1-based,
            # subtract one in order to make it correspond to the 0-based Python index
            record.seq = record.seq[sample.start - 1:sample.stop]
            # If the reverse complement is required, change the sequence accordingly
            if sample.reverse:
                record.seq = record.seq.reverse_complement()
            # Write the record to file if the file doesn't already exist
            if not os.path.isfile(sample.outputfile):
                with open(sample.outputfile, 'w') as out:
                    SeqIO.write(record, out, 'fasta')

    def __init__(self, args):
        from queue import Queue
        self.path = args.path
        self.sequencepath = os.path.join(args.sequencepath, '')
        make_path(self.sequencepath)
        self.genbankpath = os.path.join(self.path, 'genbank')
        make_path(self.genbankpath)
        self.file = os.path.join(self.path, args.file)
        self.start = args.start
        self.metadata = list()
        self.queue = Queue(maxsize=5)
        self.excelparse()

if __name__ == '__main__':
    # Argument parser for user-inputted values, and a nifty help menu
    from argparse import ArgumentParser
    import time
    # Parser for arguments
    parser = ArgumentParser(description='Download sequences from Genbank')
    parser.add_argument('path',
                        help='Specify input directory')
    parser.add_argument('-s', '--sequencepath',
                        required=True,
                        help='Path to store downloaded sequence files')
    parser.add_argument('-f', '--file',
                        help='Name of file with: "Class", "Gene", "Accession", "Start", "Stop", "Reverse" as'
                             'the headers')
    # Get the arguments into an object
    arguments = parser.parse_args()

    # Define the start time
    arguments.start = time.time()
    # Run it
    Download(arguments)
    # Print a bold, green exit statement
    print('\033[92m' + '\033[1m' + "\nElapsed Time: %0.2f seconds" % (time.time() - arguments.start) + '\033[0m')
