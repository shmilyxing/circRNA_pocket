INTRODUCTION

To use this script,perl and R should be installed. R(https://www.r-project.org/)

**for_pos

*for_pos/query_pos.pl
The Perl script for_pos/query_pos.pl is used to show the gene struct on a query site in two format,circle and linear.
You can run the script by execute comand:
perl ../query_pos.pl

*for_pos/pos_config.txt
Add your query site into for_pos/pos_config.txt.
The corresponding genome annotation file is also needed and the file path should be added into for_pos/pos_config.txt.

*for_pos/pos_format.txt
for_pos/pos_format.txt is to control the traits of output picture,such as color and width.



**for_seq

*for_seq/query_seq.pl
To use for_seq/query_seq.pl,blast should be installed.BLAST(ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+)
The Perl script for_seq/query_seq.pl is used to determine whether the PCR product Sanger sequences and RNA-Seq reads support the back splicing event.
You can run the script by execute comand:
perl ../query_seq.pl

*for_seq/config.txt
Add the file path of query sequences into for_seq/config.txt.
The genome reference and corresponding genome annotation file is also needed and their path should be added into for_seq/config.txt.The blast index of genome should be made by makeblastdb beforehand.

*for_seq/format.txt
for_seq/format.txt is to control the traits of output picture,such as color and width.



**for_sponge

*for_sponge/query_sponge.pl
To use for_sponge/query_sponge.pl,RNAhybrid or miranda should be installed.
RNAhybrid(https://bibiserv.cebitec.uni-bielefeld.de/rnahybrid) miranda(http://www.microrna.org/)
The function is not suitable for plants.
The Perl script for_sponge/query_sponge.pl is used to microRNA target on the gene struct of a query site.
You can run the script by execute comand:
perl ../query_sponge.pl

*for_sponge/sponge_config.txt
Add your query site into for_sponge/sponge_config.txt.
The corresponding mature microRNA sequences file is needed and the file path should be added into for_sponge/sponge_config.txt
You can get the mature microRNA sequences of the species you queried from miRBase(http://www.mirbase.org/cgi-bin/browse.pl?org=hsa)
The genome reference and corresponding genome annotation file is also needed and their path should be added into for_sponge/sponge_config.txt.The blast index of genome should be made by makeblastdb beforehand.

*for_sponge/sponge_format.txt
for_sponge/sponge_format.txt is to control the traits of output picture,such as color and width.


In each folder,the output example is under output_example.
The example genome reference in for_pos and for_seq is rice and in for_sponge is human(hg38).

If you upload the query sequence file to the linux server from windows operating system and some error occurred,you may use 'dos2unix' to transform