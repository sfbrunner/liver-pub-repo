import argparse
import logging
import pandas as pd
import numpy as np
import pysam
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from pile import pysam_handler
from pile import get_logger

parser = argparse.ArgumentParser(description='Provides pileup of INDELs specified in a CSV file.')
parser.add_argument('-b','--bam_file',help='Path to BAM file',required=True)
parser.add_argument('-r','--indel_csv',help='Path to CSV file containing INDELs. Columns: chrom, pos, ref, alt, type (one of Ins, Del, Complex)',required=True)
parser.add_argument('-o','--output_file',help='Path to write pileup',required=True)
#parser.add_argument('-q','--min_mapq',help='Minimal MAPQ',required=False, default=1)
#parser.add_argument('-d','--max_depth',help='Maximum depth',required=False, default=8000)
parser.add_argument('-l','--log_file',help='Path to log file',required=False, default='pileup.log')


class PileIndel(pysam_handler):

    def __init__(self, path_bam, genome_path=('/nfs/cancer_ref02/human/1000Genomes_hs37d5/genome.fa')):
        '''
        Instantiates parent class pysam_handler and establishes alignment connection.
        '''
        pysam_handler.__init__(self, genome_path=('/nfs/cancer_ref02/human/1000Genomes_hs37d5/genome.fa'))
        self.path_bam = path_bam
        self.alignment_connxn = self.connect_alignment_file(path_bam)

    def get_pileup(self, chrom, pos, ref, alt, event_type):
        '''
        Handles pileup of deletions and insertions by calling separate functions.
        '''
        #print event_type
        if event_type == 'Del':
            return self.pile_deletion(chrom, pos, ref, alt)
        elif event_type == 'Ins':
            return self.pile_insertion(chrom, pos, ref, alt)
        elif event_type == 'Complex':
            return {} #self.pile_complex(chrom, pos, ref, alt)
        else:
            return {}

    def pile_deletion(self, chrom, pos, ref, alt):
        '''
        Piles up deletions
        '''
        # Initialise
        del_seq = ref[1:]
        del_len = len(del_seq)
        pos_pre_base = pos-1
        pos_post_base = pos+del_len
        deletion_dict = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt':alt, 'depth': 0, 'type': 'Del', 'count': 0}

        # Construct the "expected" sequence, ie. Ref-DelSeq-Ref
        ref_base = str(self.fetch_seq(chromosome=chrom, start=pos, end=pos, strand='+'))
        post_base = str(self.fetch_seq(chromosome=chrom, start=pos+del_len+1, end=pos+del_len+1, strand='+'))
        expected = ref_base + del_len*'D' + post_base

        # Get pileup
        plp = self.alignment_connxn.pileup(reference = chrom, start=pos, end=pos+1, truncate=True, stepper='all')

        # Now get the observed sequences in each read
        read_dict = {}

        # Loop through reads
        for column in plp:
            for read in column.pileups:
                # Get read name
                read_name = read.alignment.query_name

                # Get observed sequence
                observed = self.get_observed_deletion(read, pos_pre_base, pos_post_base)

                # Assign the observed sequence
                read_dict[read_name] = observed

        # Evaluate whether expected and observed are identical, if so +1
        for key in read_dict:
            if read_dict[key] == expected:
                deletion_dict['count'] += 1
        
        # Calculate depth
        deletion_dict['depth'] = len(read_dict)

        # Return
        return deletion_dict

    def pile_insertion(self, chrom, pos, ref, alt):
        '''
        Piles up insertions
        '''
        # Initialise
        ins_seq = alt[1:]
        ins_len = len(ins_seq)
        pos_pre_base = pos
        pos_post_base = pos+1
        ins_dict = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt':alt, 'depth': 0, 'type': 'Ins', 'count': 0}

        # Construct the "expected" sequence, ie. Ref-InsSeq-Ref
        ref_base = str(self.fetch_seq(chromosome=chrom, start=pos_pre_base, end=pos_pre_base, strand='+'))
        post_base = str(self.fetch_seq(chromosome=chrom, start=pos_post_base, end=pos_post_base, strand='+'))
        expected = ref_base + ins_seq + post_base

        # Get pileup
        plp = self.alignment_connxn.pileup(reference = chrom, start=pos, end=pos+1, truncate=True, stepper='all')

        # Now get the observed sequences in each read
        read_dict = {}

        # Loop through reads
        for column in plp:
            for read in column.pileups:
                # Get read name
                read_name = read.alignment.query_name

                # Get observed insertion

                observed = self.get_observed_insertion(read, pos_pre_base, pos_post_base)

                # Assign the observed sequence
                read_dict[read_name] = observed

        # Evaluate whether expected and observed are identical, if so +1
        for key in read_dict:
            if read_dict[key] == expected:
                ins_dict['count'] += 1

        # Calculate depth
        ins_dict['depth'] = len(read_dict)

        # Return
        return ins_dict

    def pile_complex(self, chrom, pos, ref, alt):
        '''
        DOES NOT YET WORK!
        (Probably will never work, because Complex events require read realignment.)
        '''
        # Initialise
        #if pos==27628923:
        #    from IPython.core.debugger import Tracer; Tracer()()
        ins_seq = alt[1:]
        ins_len = len(ins_seq)
        pos_pre_base = pos
        pos_post_base = pos+1
        event_dict = {'chrom': chrom, 'pos': pos, 'ref': ref, 'alt':alt, 'depth': 0, 'type': 'Complex', 'count': 0}

        # Construct the "expected" ref sequence, ie. Ref-DelSeq-Ref
        ref_base = str(self.fetch_seq(chromosome=chrom, start=pos_pre_base, end=pos_pre_base, strand='+'))
        post_base = str(self.fetch_seq(chromosome=chrom, start=pos_post_base, end=pos_post_base, strand='+'))
        expected_ref = ref + post_base

        # Construct the "expected" alt sequence, ie. RefBase-InsSeq-RefBase
        expected_alt = alt + post_base

        # Get pileup
        plp = self.alignment_connxn.pileup(reference = chrom, start=pos, end=pos+1, truncate=True, stepper='all')

        # Now get the observed sequences in each read
        read_dict = {}

        # Loop through reads
        for column in plp:
            for read in column.pileups:
                read_name = read.alignment.query_name
                #if read_name == 'HX4_22523:3:2109:17868:3858':
                #    from IPython.core.debugger import Tracer; Tracer()()
                # Get the observed deletion event
                observed_del = self.get_observed_deletion(read, pos_pre_base, pos_post_base)

                # Get the observed deletion event
                observed_ins = self.get_observed_insertion(read, pos_pre_base, pos_post_base)

                # Assign the observed sequence
                read_dict[read_name] = {'del':observed_del, 'ins':observed_ins}

        # Evaluate whether expected and observed are identical, if so +1
        for key in read_dict:
            if read_dict[key] == expected:
                event_dict['count'] += 1

        # Calculate depth
        event_dict['depth'] = len(read_dict)

        # Return
        return event_dict

    def get_observed_deletion(self, read, pos_pre_base, pos_post_base):# Loop through all relevant query positions
        aligned_pairs = read.alignment.aligned_pairs
        relevant_queries = [pair[0] for pair in aligned_pairs if (pair[1]>=pos_pre_base and pair[1]<=pos_post_base)]

        # Get the sequences at those positions
        observed = ''
        for query in relevant_queries:
            if query == None:
                observed += 'D'
            else:
                observed += read.alignment.query_sequence[query]
        return observed

    def get_observed_insertion(self, read, pos_pre_base, pos_post_base):
        aligned_pairs = read.alignment.aligned_pairs
        first_q_base = [i for i in range(0, len(aligned_pairs)) if aligned_pairs[i][1]==(pos_pre_base-1)]
        last_q_base  = [i for i in range(0, len(aligned_pairs)) if aligned_pairs[i][1]==pos_post_base]
        if first_q_base and last_q_base:
            observed = read.alignment.query_sequence[first_q_base[0]:last_q_base[0]]
        else:
            observed = ''
        return observed


def run_pileup(path_indel, path_bam, outpath, tbl_delim=',', logger_obj=None):
    '''
    Runs pileup on positions specified in BED file, on the given BAM file.
    
    Args:
    path_csv: path to CSV file specifying INDELs (see args)
    path_bam: path to a BAM file
    outpath: path to write pileup
    tbl_delim: delimiter used in CSV file (should be ',')
    '''
    
    # Import CSV file to Pandas DF
    indel_tbl = pd.read_table(path_indel, delimiter=tbl_delim) #, escapechar='#', header=True)
    #logger_obj.info(indel_tbl.head())

    # Connect BAM file
    pile_indel = PileIndel(path_bam = path_bam)
    
    dict_lst = []
    counter = 0
    first_write = 1

    # Loop through CSV coordinates
    for index, row in indel_tbl.iterrows():
        counter += 1
        if logger_obj:
            if (counter % 100) == 0:
                logger_obj.info('Piling INDEL {c} of {t}'.format(c=str(counter), t=indel_tbl.shape[0]))
        row_dict = dict()
        
        try:
            chrom = str(row['chrom'])
            pos = row['pos']
            ref = row['ref']
            alt = row['alt']
            event_type = row['type']
            plp = pile_indel.get_pileup(chrom, pos, ref, alt, event_type)
        except Exception, e:
            print('Pileup skipped with error message: {err}'.format(err=str(e)))
        else:
            if plp:
                this_tbl = pd.DataFrame(plp, index=[0])
                #from IPython.core.debugger import Tracer; Tracer()()
                this_tbl = this_tbl[['chrom', 'pos', 'ref', 'alt', 'type', 'depth', 'count']]
                if first_write==1:
                    this_tbl.to_csv(outpath, index=False)
                    first_write=0
                else:
                    this_tbl.to_csv(outpath, index=False, mode='a', header=False)
    
    return True

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
    run_pileup(path_indel = args.indel_csv, path_bam = args.bam_file, outpath=args.output_file, logger_obj=logger)
    logger.info('Pileup complete, storing output to disk at {p}'.format(p=args.output_file))