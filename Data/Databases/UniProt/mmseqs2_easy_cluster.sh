mkdir tmp

PREFIX="Arabidopsis_thaliana_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Caenorhabditis_elegans_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Danio_rerio_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Drosophila_melanogaster_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Homo_sapiens_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Mus_musculus_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Saccharomyces_cerevisiae_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

PREFIX="Xenopus_tropicalis_proteome"
mmseqs easy-cluster $PREFIX.fa.gz $PREFIX.nr20 tmp --min-seq-id 0.2 -s 7.5
gzip -c $PREFIX.nr20_rep_seq.fasta > $PREFIX.nr20.fa.gz
rm $PREFIX.nr20_*

rm -rf tmp