cut -d "," -f 4 Arabidopsis_thaliana.csv | cut -d "\"" -f 2 > Arabidopsis_thaliana_TFs.txt
cut -d "," -f 4 Caenorhabditis_elegans.csv | cut -d "\"" -f 2 > Caenorhabditis_elegans_TFs.txt
cut -d "," -f 4 Drosophila_melanogaster.csv | cut -d "\"" -f 2 > Drosophila_melanogaster_TFs.txt
cut -d "," -f 4 Homo_sapiens.csv | cut -d "\"" -f 2 > Homo_sapiens_TFs.txt
cut -d "," -f 4 Mus_musculus.csv | cut -d "\"" -f 2 > Mus_musculus_TFs.txt
cut -d "," -f 4 Saccharomyces_cerevisiae.csv | cut -d "\"" -f 2 > Saccharomyces_cerevisiae_TFs.txt
