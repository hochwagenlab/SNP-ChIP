# README
# ChIP-seq Pipeline version 3 
( updated: 1/22/18 )

## UPDATES
1. This update is to allow for paired-end sequencing analysis.
2. Two new scripts have been added to the folder: Bowtie_PE_B4.sbatch
   and MACS2_FE_W4.sbatch. They are not currently called by 
   ChIPseq-pipeline_v3.sbatch.
3. See ChIPseq_PE_example.sh for a template to run the paired-end
   sequencing pipeline.


## PREPARATION
1. Move entire folder into your home on prince (keep folder name as
   "ChIPseq_Pipeline_v3")
2. Move the 4 fasta files into a folder called "Library" in your home on prince
3. Open each [4] sbatch script and adjust the email address
    * You can do it from the CLI by running the following command inside
the directory (change "X" to your user name; this will replace "tem298"):

    `find . -type f | xargs sed -i 's/tem298@nyu.edu/X@nyu.edu/g'`

4. Start running jobs (sbatch ~/ChIPseq_Pipeline/ChIPseq-pipeline_v3.sbatch)
5. Keep an eye out for the following errors:  
   **slurmstepd: error: Exceeded step memory at some point.**  
   - If you see this in association with a Bowtie job:  
     (should not be an issue if files are less than 1.2Gb in zipped format)  
     - adjust the memory requirements for Bowtie_B3.sbatch   
     - delete all associated fastq, fastq.gz, and sam files   
     - rerun job  
   - If you see this error with a MACS2 job:  
     - ensure all .bdg and .wig files are zipped    
     - if not, rerun job with more memory for MACS2_FE_W3.sbatch

   **slurmstepd: error: Exceeded job memory at some point.**  
   - If you see this error:  
     - delete all associated files   
     - increase memory on associated sbatch job (Bowtie or MACS2)   
     - rerun job

## SOME RUN EXAMPLES
Example 1: Start with input file A.fastq and ChIP file B.fastq, map to SacCer3, 
and create bedgraph files and both narrow and broad peak files 
(will also create wiggle plot of the rDNA)
```Bash
sbatch --export INPUT=A.fastq,CHIP=B.fastq,GEN="SacCer3",TAGI="A",TAGC="B",PEAK="BOTH" \ 
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

Example 2: Get broad peak files for replicate input files A_1.sam and A_2.sam, and
ChIP files B_1.sam and B_2.sam [Note change from previous pipeline.]
```Bash
sbatch --export INPUT="A_1.sam A_2.sam",CHIP="B_1.sam B_2.sam",REP="B",FLMKR="1-3",PEAK="BROAD",WIG="F" \
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

Example 3: Get normalized wiggle files (and narrow peaks) with ChIP file A1.sam, input file A2.sam,
mock ChIP file B1.sam, and mock input file B2.sam
```Bash
sbatch --export INPUT=A2.sam,CHIP=A1.sam,INCON=B2.sam,CHCON=B1.sam,FLMKR="1-2" \
~/ChIPseq_Pipeline_v3/ChIPseq-pipeline_v3.sbatch
```

## OUTPUT FILE NAMES
All outputs will get a version name identifier:
  - B3 means that the Sam file was created with this pipeline version
  - W3 means that the peak/bedgraph/wiggle files were made with this pipeline version
  - B3W3 means that the everything was made with this pipeline version

### SAM:
All new sam files will be given the following names: $TAG-$GENROOT.sam
where $TAG is the user defined name and $GENROOT is one of the following:
  - SK1Yue-PM_B3
  - SK1Yue-2mis_B3
  - SacCer3-2mis_B3
  - SK1K-PM_B3
  - SacCer3-rDNA_B3
  - SK1K-rDNA_B3

### MACS2 FOLDER:
Macs2 files will now go in folders with the notation: $TAGC-$GEN-$VER-MACS2 where  
  $TAGC is either the tagname of the ChIP file or the user-defined name $REP,  
  $GEN is either SacCer3 or SK1K, and  
  $VER will currently be either B3W3, 2mis_W3, or PM_W3.

$FLMKR will be incorporated into all MACS2 output files.  
The term "Reps" will be incorporated into all MACS2 output files when applicable.

### MACS2 FILES:
See MACS2 FOLDER with the following changes:
  - All files will include the method of Bowtie analysis: 2mis, PM, or rDNA.
  - Wiggle and bedgraph files will include the normalization method: FE.
  - ChIP vs ChIP analysis is identified as either: ChvCh or InvIn.

## ERROR/OUTPUT FILES
- The pipeline automatically makes two error/output files that I haven't figured  
  out how to remove yet. Both should be essentially empty at all times. One   
  (closing_ChIPseq_$JOBID.out) will go to /scratch/$USER and the other (slurm_$JOBID.out)  
  will go to the working directory where you submit the job.  
- The pipeline will initially make a /scratch/$USER/MACS2_pipeline_$CHIP.txt.  
  If none of the bowtie or macs2 jobs are initiated, this file will include  
  all necessary error messages.
- While running, the pipeline will create many .out files in /scratch/$USER. These   
  will eventually be consolidated along with the above .txt into a single file called  
  one of the following based upon the given arguments by order of preference:  
       - $REP_$JOBID.out
       - $TAGC_$JOBID.out
       - $CHIP_$JOBID.out

## FOR TROUBLESHOOTING
On the off chance a set of files does not get consolidated into a single file:
1) Check the file /scratch/$USER/MACS2_pipeline_$CHIP.txt. It will tell   
   you all of the jobs that ran for this pipeline run as a colon-separated list.
2) For each job, just type the following command to see what might have gone wrong: 
```Bash
less /scratch/$USER/*[insert JobID here]* 
```  
 These may be called:  
   - Bowtie_$JOBID.out  
   - MACS2_FE_$JOBID.out  
   - closing_ChIPseq_$JOBID.out  

## NEW FEATURES
- Bowtie_B3.sbatch automatically checks for a file in both zipped and unzipped   
  forms, and proceeds accordingly. Therefore, inputs into the pipeline can be   
  written in either format.
- Bowtie_B3.sbatch includes a check for the databases needed, and if not found,   
  attempts to build them. No need for a separate bowtie-build.
- When starting with two (or four for ChIP vs ChIP normalization) fastq files,   
  the pipeline automatically processes the rDNA along with the rest of the pipeline.
- For details on the Bowtie or MACS2 conditions, see Bowtie_B3.sbatch or  
  MACS2_FE_W3.sbatch.

## KNOWN ISSUES
- There are currently no error checks for the following:
  - Do all SAM files have the same root?
  - Have wiggle or peak files already been created?
  - Is the FileMaker ID currently in the name? (Currently there is only a   
    superficial search).
  - Bowtie_B3.sbatch still takes 5 minutes even if it should only  
    be checking file names and unzipping/zipping a single fastq file.
  - Two extra slurm output files are created
  - No simple method to create an rDNA wiggle file from two  
       (or more) sam files.
