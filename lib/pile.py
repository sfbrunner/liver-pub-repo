import argparse
import logging
import pandas as pd
import numpy as np
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

parser = argparse.ArgumentParser(description='Creates pileup tables from BED file locations.')
parser.add_argument('-b','--bam_file',help='Path to BAM file',required=True)
parser.add_argument('-r','--bed_regions',help='Path to BED file',required=True)
parser.add_argument('-o','--output_file',help='Path to write pileup',required=True)
parser.add_argument('-q','--min_mapq',help='Minimal MAPQ',required=False, default=1)
parser.add_argument('-d','--max_depth',help='Maximum depth',required=False, default=8000)
parser.add_argument('-l','--log_file',help='Path to log file',required=False, default='pileup.log')


class pysam_handler():
    '''
    Provides basic functionality to interact with the pysam package.
    '''

    def __init__(self, genome_path=('/nfs/cancer_ref02/human/1000Genomes_hs37d5/genome.fa')):
        '''
        Args:
        - genome_path: path to FASTA file containing reference genome sequence.
                        default: ('/nfs/cancer_ref02/human/1000Genomes_hs37d5/genome.fa').
                        An index (.fai) needs to exist in the same directory.
        '''
        self.genome_path = genome_path
        self.genome = pysam.FastaFile(self.genome_path)

        # Initialise helper variables
        self.strand_allowed_vals = ['-', '+']

    def fetch_seq(self, chromosome, start, end, strand = '+'):
        '''
        Fetch DNA sequence from genome file

        Args:
        - chromosome: chromosome, needs to be compatible with reference names
                        contained in self.genome
        - start: start of sequence to fetch (1-based)
        - end: end of sequence to fetch (1-based)
        - strand: strand of sequence to fetch, takes '+' or '-' (default: '+')
        '''
        if strand not in self.strand_allowed_vals:
            return ''
        elif strand == '+':
            return Seq(self.genome.fetch(reference=chromosome, start=start-1, end=end))
        elif strand == '-':
            return Seq(self.genome.fetch(reference=chromosome, start=start-1, end=end)).reverse_complement()

    def connect_alignment_file(self, fpath_alignment_file):
        '''
        Returns a connection to an alignment file (BAM/CRAM/SAM)

        Args:
        fpath_alignmet_file: file path to alignment file
        '''
        return pysam.AlignmentFile(fpath_alignment_file, 'r'+self.check_alignment_read_mode(fpath_alignment_file))

    def check_alignment_read_mode(self, fpath):
        '''
        Determines whether alignment file is CRAM, BAM or SAM and sets read mode
        relevant to pysam accordingly

        Args:
        fpath: path to bam file
        '''

        # Determine required read mode based on file extension
        file_ext = fpath.split('.')[-1]
        if file_ext == 'bam':
            read_mode = 'b'
        elif file_ext == 'cram':
            read_mode = 'c'
        else:
            read_mode = ''

        return read_mode

    def get_pileup_from_region(self, alignment_connxn, chrom, start_pos, end_pos, max_depth=None, min_mapq=-1):
        '''
        Returns a Pandas table containing a pileup over the defined region.

        Args:
        alignment_connxn: connection to alignment file, as returned by method
            connect_alignment_file()
        chrom: chromosome of region of interest
        start: start of region of interest
        end: end of region of interest
        max_depth: maximal depth permitted at position. If depth is greater, then ValueError is thrown.
        min_mapq: any read with mapq > min_mapq will be used
        '''

        # Initialise list of dictionaries
        lst_of_pos_dict = []

        # Initialise pileup
        #print 'Chr: ' + str(chrom) + 'Start position: '+str(start_pos)+', end position: '+str(end_pos)
        #plp = alignment_connxn.pileup(reference = chrom, start=start_pos, end=end_pos, truncate=True) # Do not use stepper='samtools', that would exclude any reads mapped in unexpected orientations. Those can be interesting, although may at times represent noise.
        plp = alignment_connxn.pileup(reference = chrom, start=start_pos, end=end_pos, truncate=True, stepper='all')

        # Loop through positions
        for column in plp:
            # Extract position and reference base
            ref_pos = column.pos
            ref_base = str(self.fetch_seq(chromosome=chrom, start=ref_pos, end=ref_pos, strand='+'))
            mapq_lst = []
            
            # Initialise position dictionary
            pos_dict = {'chrom':chrom, 'pos': ref_pos, 'ref': ref_base, 'A':0, 'C':0, 'G':0, 'T':0, 'N':0, 'INS':0, 'DEL':0, 'a':0, 'c':0, 'g':0, 't': 0, 'n':0, 'ins':0, 'del':0, 'mapq':0}
            
            # Loop through reads
            counter=0
            for read in column.pileups:
                # Discard if soft-clipped
                if read.alignment.get_reference_positions(full_length=True)[read.query_position-1]:
                    
                    # Discard if mapq below cutoff
                    this_mapq = read.alignment.mapping_quality
                    if this_mapq<min_mapq:
                        continue
                    mapq_lst.append(read.alignment.mapping_quality)

                    # Check if indel
                    if read.indel > 0:
                        query_base = 'INS'
                        #from IPython.core.debugger import Tracer; Tracer()()
                    elif read.is_del == 1:
                        query_base = 'DEL'
                    else:
                        query_base = read.alignment.query_sequence[read.query_position-1] # Pysam is 0-based
                    
                    # Check strand
                    if read.alignment.is_reverse:
                        query_base = query_base.upper()
                    else:
                        query_base = query_base.lower()

                    # Augment base count
                    pos_dict[query_base] += 1
                    counter += 1

                    # DEBUGGING
                    #print('Read {r:s}, base: {b:s}'.format(r=read.alignment.qname, b=query_base))
                
                    if max_depth:
                        if counter > max_depth:
                            raise ValueError('Reached maximum depth at position {c}:{p}'.format(c=str(chrom), p=str(ref_pos)))
            if mapq_lst:
                pos_dict['mapq'] = np.median(np.array(mapq_lst))
            else:
                pos_dict['mapq'] = 0
            lst_of_pos_dict.append(pos_dict)

        # Create Pandas table from list of dictionaries
        df = pd.DataFrame(lst_of_pos_dict)
        #print df
        if not df.empty:
            df = df[['chrom', 'pos', 'ref', 'A', 'C', 'G', 'T', 'N', 'INS', 'DEL', 'a', 'c', 'g', 't', 'n', 'ins', 'del', 'mapq']] # Reorder columns
        return df

def run_pileup(path_bed, path_bam, outpath, max_depth = 8000, min_mapq=1, bed_delim='\t', logger_obj=None):
    '''
    Runs pileup on positions specified in BED file, on the given BAM file.
    
    Args:
    path_bed: path to a BED file
    path_bam: path to a BAM file
    max_depth: maximal depth to run pileup
    min_mapq: minimal mapq required to count a base
    bed_delim: delimiter used in BED file
    '''
    
    # Import BED file to Pandas DF
    bed_regions = pd.read_table(path_bed, delimiter=bed_delim, escapechar='#', header=None)
    bed_regions.columns = ['chr', 'start', 'end']
    #bed_regions = bed_regions.loc[1:1000,:]
    
    # Connect BAM file
    ph = pysam_handler()
    bam_connxn = ph.connect_alignment_file(path_bam)
    
    dict_lst = []
    counter = 0
    first_write = 1

    # Loop through BED coordinates
    for index, row in bed_regions.iterrows():
        counter += 1
        if logger_obj:
            if (counter % 100) == 0:
                logger_obj.info('Piling BED row {c} of {t}'.format(c=str(counter), t=bed_regions.shape[0]))
        row_dict = dict()
        try:
            this_chrom = str(row['chr'])
            this_start = row['start']
            this_end   = row['end']
            plp = ph.get_pileup_from_region(bam_connxn, this_chrom, this_start, this_end+1, min_mapq=min_mapq)
        except Exception as e:
            print('Pileup ({c}:{p}) skipped with error message: {err}'.format(c=str(this_chrom), 
                                                                                                p=str(this_start),
                                                                                                err=str(e)))
        else:
            if not plp.empty:
                if first_write==1:
                    plp.to_csv(outpath, index=False)
                    first_write=0
                else:
                    plp.to_csv(outpath, index=False, mode='a', header=False)
    
    return True

def get_logger(log_file_path):
    logger = logging.getLogger(__name__)
    logger.setLevel(logging.DEBUG)
    # create file handler which logs even debug messages
    fh = logging.FileHandler(log_file_path)
    fh.setLevel(logging.DEBUG)
    # create console handler with a higher log level
    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    # create formatter and add it to the handlers
    formatter = logging.Formatter(fmt='%(asctime)s [%(levelname)s]: %(message)s', datefmt='%m/%d/%Y %I:%M:%S %p')
    fh.setFormatter(formatter)
    ch.setFormatter(formatter)
    # add the handlers to the logger
    logger.addHandler(fh)
    logger.addHandler(ch)

    return logger

if __name__ == '__main__':
    args = parser.parse_args()
    
    # Initialise logger
    logger = get_logger(log_file_path = args.log_file)
            
    # Start logging
    logger.info('\n------\nSetting up pileup job \n------\n')
    logger.info('Created log at path {path}.'.format(path=args.log_file))
    logger.info('User-defined arguments:')
    for arg, value in sorted(vars(args).items()):
        logger.info("Argument %s: %r", arg, value)
    
    # Run the job
    logger.info('Running pileup.')
    run_pileup(path_bed = args.bed_regions, path_bam = args.bam_file, outpath=args.output_file,
               max_depth = args.max_depth, min_mapq=args.min_mapq, logger_obj=logger)
    logger.info('Pileup complete, storing output to disk at {p}'.format(p=args.output_file))
    #plp_tbl.to_csv(args.output_file, index=False)


#path_bam = "/nfs/users/nfs_s/sb50/data/1666/cancer_ref01/PD36715b_lo001/PD36715b_lo001.sample.dupmarked.bam"
#path_bed = "~/mylustre2/1666/shearwater_long/unique_pos_vaf15.bed"

#run_pileup(path_bed, path_bam, bed_delim=' ')
