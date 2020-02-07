zgrep "^AC   " *_TFs.txt.gz | head | perl -e '
%uniprot;
while(<>){
    chomp;
    @matches=($_=~m/(\S+);/g);
    foreach $match(@matches){
        $uniprot{$match}=1;
    }
}
# foreach $uniacc(keys %uniprot){
#     print "$uniacc\n";
# }
open(FH, "<", "../Pfam/Pfam-A.regions.uniprot.tsv");
while(<FH>){
    ($uniacc)=$_=~m/^(\S+)/;
    if(defined $uniprot{$uniacc}){
        print $_;
    }
}
close(FH);
'