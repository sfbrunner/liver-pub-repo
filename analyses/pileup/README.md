# Running pileups using pilepy (pysam)

## Install Python 2.7 virtualenv
* Install virtualenv using instructions specific to your OS. Within the Sanger/CASM infrastructure, follow the steps on Confluence: https://confluence.sanger.ac.uk/display/IT/Python+virtualenv+setup+guide
* Activate the virtualenv and check that pip is the virtualenv's pip (which pip)
* Install the necessary packages using requirements.txt. Command: `pip install -r requirements.txt`

## Run SNV pileup
Activate your virtualenv and use the following code example to run an SNV pileup:
```
python ../../lib/pile.py \
    --bam_file /nfs/cancer_ref01/nst_links/live/1666/PD37105b_lo001/PD37105b_lo001.sample.dupmarked.bam \
    --bed_regions input/example.bed \
    --output_file output/pilepy_example.csv \
    --log_file output/pilepy_example.log
```

### Run SNV pileup using FARM cluster
A simple bsub script is provided in `snv_pileup.sh`. The script reads the locations of BAM files from a text file (here given in `input/all_bams.txt`) and runs an array job over all the BAM files, ie. running one parallel task per BAM file. To run your own array job, copy and modify the `snv_pileup.sh` file, then submit with `bsub < [your_file].sh`. 