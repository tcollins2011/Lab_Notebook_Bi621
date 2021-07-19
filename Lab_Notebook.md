<mark>General Notes</mark>

-------
2021/7/8
-------
<!-- Way to set up interactive node -->
srun --account=bgmp --partition=bgmp --nodes=1 --ntasks-per-node=1 --time=2:00:00 --cpus-per-task=1 --pty bash

<!-- Reminder of useful conda commands to set up environment -->
conda evn list
conda activate bgmp_py39


<mark>PS6</mark>

------------------
7/14/21
------------------

<!-- Velvet Reports things in terms of k-mers. So to convert back to read coveage you need to un k-merize it.  -->
contig length = kmer  + reported length - 1
read coverage = kmer_coverage / kmer


<!-- It turns out that weighted averages are important and that I should be using them on assembly data-->
sum(C * L) / sum(L)

-------------------
7/15/21
-------------------

<!-- Fun error with velvet and kmer lenth 41.  -->

    I have removed my colorful language, but there still seems to be an issue with using the separate tag on paired end reads. According to the documentation this should be correct but it is still occasionally resulting in contigs on k41 that are billions of bp long. This could require additional research as even adding Pete's command hasn't seem to totally fix it.

-------------------
7/16/21
-------------------
<!-- Run time max information for Velvet Bash script -->

<!-- Velvet H -->
	User time (seconds): 49.92
	System time (seconds): 1.05
	Percent of CPU this job got: 337%
	Maximum resident set size (kbytes): 1099756

 <!-- Velvet G -->
	User time (seconds): 96.15
	System time (seconds): 0.85
	Percent of CPU this job got: 307%
	Maximum resident set size (kbytes): 569468

    This should be a good estimate for how many cpus and how much ram I need to rerun my assembly sequences.

-------------------
7/18/21
-------------------
<!-- Fixed the issue that was present in matplot lib when looping over multiple graphs -->
    add the below commands to fully reset matplot lib inbetween runs so that the graphs don't hold previous info
        plt.clf()
        plt.cla()

<!-- Changed the program structure. -->
    ps6.py will now run the unit test when run by default

    main.py will run analysis on all of the target contig files when called now. It will also create all necessary graphs
    and tables as well.

<!-- How to run from scratch -->
    Run the ps6.py file to test that it works
    Run the run_velvet.sh script to generate the 7 different runs 
    Run main.py to generate a csv and all graphs, and tsv files
    If you want a nice markdown table use a convertor online for the csv table


<mark>PS7</mark>

---------------
2021/7/14
---------------
I'm over complicating my fasta function for no reason but it is finally passing its tests and working so I will leave it. 

The whole code works but I'm am getting 23,478 proteins for humans and not the correct 23,477

 <!-- This is now solved. I just needed to use a dictionary key value reassignment instead of an update  -->


---------------
2021/7/15
---------------
The whole code is now functional and returning the correct stats files as compared with my peers. At this point althoguh its in a strange dictionary state I can move on to Blast P and analysis in the morning.


---------------
2021/7/16
---------------
<!-- Blast P Output Info -->

	Command being timed: "makeblastdb -in human_longest.fa -out human_blast.db -dbtype prot -parse_seqids"
	User time (seconds): 0.61
	System time (seconds): 0.01
	Percent of CPU this job got: 87%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.72
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 15432

    Command being timed: "makeblastdb -in danio_longest.fa -out danio_blast.db -dbtype prot -parse_seqids"
	User time (seconds): 0.78
	System time (seconds): 0.02
	Percent of CPU this job got: 87%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:00.92
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 16788

	Command being timed: "blastp -query human_longest.fa -evalue 1e-6 -use_sw_tback -db danio_blast.db -num_threads 8 -outfmt 6 -out homo_sapien_vs_danio_db.txt"
	User time (seconds): 22641.71
	System time (seconds): 11.91
	Percent of CPU this job got: 161%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:54:02
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 118360
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 5961577

    	Command being timed: "blastp -query danio_longest.fa -evalue 1e-6 -use_sw_tback -db human_blast.db -num_threads 8 -outfmt 6 -out danio_vs_human_db.txt"
	User time (seconds): 22928.20
	System time (seconds): 12.81
	Percent of CPU this job got: 156%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 4:03:55
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 102188
	Average resident set size (kbytes): 0
	Major (requiring I/O) page faults: 0
	Minor (reclaiming a frame) page faults: 5793185

------------
2021/7/18
------------
<!-- To run from scratch -->
run long_protein.py on the target fasta file
run blastp.sh on the newly generated fasta file
run summary.sh on a sam file 

<mark>PS8</mark>

------------
2021/7/8
------------

tarted PS8 and created its directory structure.

running conda enironment: bgmp_py39 

Make sure to switch to the interactive node and not stay on the login node. I keep checking the node after I run something ...

Pulled two Zebra Danio files from Enbemble: 

Danio_rerio.GRCz11.dna.primary_assembly.fa.gz
Danio_rerio.GRCz11.104.gtf.gz 

<!-- Useful command for keeping track of resources used -->

/usr/bin/time

<!-- Remember to set up a sbatch script like this -->
#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8


<!-- Script used to generate star database -->

#!/bin/bash

#SBATCH --account=bgmp
#SBATCH --partition=bgmp
#SBATCH --cpus-per-task=8

conda activate bgmp_py39

/usr/bin/time -v STAR --runThreadN 8 --runMode genomeGenerate \
--genomeDir /projects/bgmp/tcollin2/bioinfo/Bi621/PS/ps8-tcollins2011/Danio_rerio.GRCz11.dna.ens104.STAR_2.7.9a \
--genomeFastaFiles /projects/bgmp/tcollin2/bioinfo/Bi621/PS/ps8-tcollins2011/dre/Danio_rerio.GRCz11.dna.primary_assembly.fa \
--sjdbGTFfile /projects/bgmp/tcollin2/bioinfo/Bi621/PS/ps8-tcollins2011/dre/Danio_rerio.GRCz11.104.gtf


 exit



 <!-- Star align reads stats -->
 Command being timed: "STAR --runThreadN 8 --runMode alignReads --outFilterMultimapNmax 3 --outSAMunmapped Within KeepPairs --alignIntronMax 1000000 --alignMatesGapMax 1000000 --readFilesCommand zcat --readFilesIn /projects/bgmp/shared/Bi621/dre_WT_ovar12_R1.qtrim.fq.gz /projects/bgmp/shared/Bi621/dre_WT_ovar12_R2.qtrim.fq.gz --genomeDir /projects/bgmp/tcollin2/bioinfo/Bi621/PS/ps8-tcollins2011/Danio_rerio.GRCz11.dna.ens104.STAR_2.7.9a --outFileNamePrefix Danio_rerio"
	User time (seconds): 1424.96
	System time (seconds): 11.55
	Percent of CPU this job got: 717%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 3:20.16
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 15886368


-----------
2021/7/12
-----------
    Inital star commands seemed to work and all that is needed is to figure out the blast documenation and to pase the file

------------
2021/7/14
------------

<!-- Stats for running blastp -->

	Command being timed: "samtools sort /home/tcollin2/bgmp/bioinfo/Bi621/PS/ps8-tcollins2011/Aligned/Danio_rerioAligned.out.bam"
	User time (seconds): 52.43
	System time (seconds): 1.03
	Percent of CPU this job got: 99%
	Elapsed (wall clock) time (h:mm:ss or m:ss): 0:53.91
	Average shared text size (kbytes): 0
	Average unshared data size (kbytes): 0
	Average stack size (kbytes): 0
	Average total size (kbytes): 0
	Maximum resident set size (kbytes): 917976

------------
2021/7/18
------------
<!-- Potential problem in python script -->
    So I was uncertain if we should discard all paired end reverse reads or only discard them if the forward reads didn't map but the reverse reads did. For the moment I'm just discarding all reserve end reads. This may need to be fixed.

<!-- How to run everything from scratch -->
    run Star_Database.sh
    run Star_align.sh
    run samtools_report.sh
    run main.py

