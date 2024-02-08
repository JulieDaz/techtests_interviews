### biographica tech test

## versions
## blastp version 2.9.0+
## python3 version 3.7.1

## create a directory to output the res
mkdir res

# ## remove the dot at the end of the sequences in assembly_1.prot.fa
sed s/.$//g data2/assembly_1.prot.fa > data2/assembly_1.prot.dotrm.fa

# ## Step1: use blast to retrieve identical/almost identical sequences between the 2 assemblies
makeblastdb -dbtype prot -in data2/assembly_1.prot.dotrm.fa 
blastp -query data2/assembly_2.prot.fa -db data2/assembly_1.prot.dotrm.fa -evalue 10e-3 -out res/blast_assembly12_output.tsv -num_threads 7 -outfmt "6 qseqid sseqid pident length mismatch gapopen qlen qstart qend slen sstart send evalue bitscore" 

# Step2: use python to parse blast output and create the output of the desambiguated sequences
python3 submission/biographica_task1.py --fasta1 data2/assembly_1.prot.fa --fasta2 data2/assembly_2.prot.fa --matrix res/blast_assembly12_output.tsv --out res/assembly1_2_desambiguate.prot.fa

## Step3: Retrieve the proteins with the activity, acting on ester bonds
python3 submission/biographica_task2.py --fasta1 data2/assembly_1.prot.fa --fasta2 data2/assembly_2.prot.fa --merged res/assembly1_2_desambiguate.prot.fa --out res/assembly1_2_proteins_hydrolase.fa

