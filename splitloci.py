#!/usr/bin/env python

import pysam
import plac
from itertools import *
import os
import os.path
import shutil
import logging

def prepare_output_dir(dir, clobber=False):
    if os.path.exists(dir):
        if clobber:
            logging.warn("Clobbering existing output directory: %s", dir)
            shutil.rmtree(dir)
        else:
            raise Exception("Output directory already exists: %s", dir)
    logging.debug("Creating dir: %s", dir)
    os.mkdir(dir)

class DummyFile(object):
    """A dummy object with a close() method."""
    def close(self):
        pass
    def write(self, *args, **kwargs):
        raise NotImplementedError("Cannot write to dummy file")

def fetch_mapped(sam, *args, **kwargs):
    for read in sam.fetch(*args, **kwargs):
        if not read.is_unmapped:
            yield read

def get_pos(read):
    return (read.tid, read.pos)

def get_locus_limit(read):
    if read.is_proper_pair and \
       not read.mate_is_unmapped and \
       read.tid == read.mrnm and \
       read.mpos > read.pos:
        pos = read.pos + read.isize
    else:
        pos = read.aend
    return (read.tid, pos)

@plac.annotations(
    # arg=(helptext, kind, abbrev, type, choices, metavar)
    sam=("Output in SAM format instead of BAM.", "flag", "S"),
    quiet=("Do not print informational messages.", "flag", "q"),
    verbose=("Print debug messages that are probably only useful if something is going wrong.", "flag", "v"),
    clobber_output_directory=("If output dir exists, delete it and all of its contents to make way for the output.", "flag", "C"),
    min_reads_per_file=("Minimum number of reads to put in each file. With this option, multiple smaller loci will be collected into a single file until that file reaches the limit.", "option", "m", int, None, 'N'),
    infile=("Input file, in either SAM or BAM format. Must be sorted by reference position and have proper headers.", "positional"),
    output_dir=("Name of a directory to create and place output files into.", "positional"),
    )
def main(sam, quiet, verbose, clobber_output_directory, infile, output_dir, min_reads_per_file=1):
    """Split a SAM file into (approximately) one file per locus.

    A locus is a contiguous region with non-zero coverage, surrounded
    by regions of zero coverage. Splitting into loci allows one to
    parallelize a Cufflinks assembly across machines in a cluster.

    The output directory will contain files named like loc000001.bam,
    loc000002.bam, and so on, each containing one full locus (or
    several full loci when using -m). Each output file will have an
    identical header to the input file."""
    # Configure logging
    if quiet:
        logging.basicConfig(level=logging.WARN)
    elif verbose:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)

    # Report options
    logging.debug("Options: %s" % {
        'sam': sam,
        'quiet': quiet,
        'verbose': verbose,
        'clobber_output_directory': clobber_output_directory,
        'infile': infile,
        'output_dir': output_dir,
        'min_reads_per_file': min_reads_per_file,
    })

    # Set up input
    sam_input = pysam.Samfile(infile)
    reads = fetch_mapped(sam_input)

    # Set up output
    prepare_output_dir(output_dir, clobber_output_directory)
    locfile_template = os.path.join(output_dir, "loc%06d")
    if sam:
        locfile_template += ".sam"
        output_mode = "wh"
    else:
        locfile_template += ".bam"
        output_mode = "wb"
    def locfile(n):
        """Return the nth locfile"""
        return pysam.Samfile(locfile_template % n, output_mode, template=sam_input)

    # Initialize loop variables
    loc_count = 1
    logging.info("Starting locus file %s" % loc_count)
    output = locfile(loc_count)
    # Need the first read in order to initialize loc_limit
    first_read = next(reads)
    output.write(first_read)
    loc_limit = get_locus_limit(first_read)
    read_count = 1

    for read in reads:
        assert read.tid >= loc_limit[0]
        if (get_pos(read) > loc_limit):
            if read_count < min_reads_per_file:
                logging.info("Adding another locus to locus file %s" % loc_count)
            else:
                # Need to start a new locus
                logging.info("Finished locus file %s with %s reads" % (loc_count, read_count))
                output.close()
                loc_count += 1
                logging.info("Starting locus file %s" % loc_count)
                output = locfile(loc_count)
                # Reset read count for the new file
                read_count = 0
        output.write(read)
        read_count += 1
        # Extend loc_limit to the end of the fragment represented by the read
        loc_limit = max(loc_limit, get_locus_limit(read))
    logging.info("Finished locus file %s with %s reads" % (loc_count, read_count))
    output.close()

# Entry point
def plac_call_main():
    return plac.call(main)

if __name__=="__main__":
    plac_call_main()
