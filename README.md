# barcode-constrained-phylogeny
This repository contains code and data for building topologically-constrained phylogenies.

NEEDS TRANSLATING! Also need to add the two files to repo?

De internationale database BOLD heeft DNA barcodes voor honderdduizenden soorten (met vaak meerdere barcodes per soort). 
Als je de data in BOLD zou filteren op merker, er zitten immers meerdere merkers in de database, en die vervolgens per merker zou alignen, dan zou dit redelijke data zijn om fylogenetische bomen mee te bouwen. Echter, er zijn twee beperkingen: 1) er zit een bovengrens aan de grootte van alignments en bomen die te bouwen zijn, en 2) de barcodes zijn vrij kort dus geven vrij weinig signaal om grote bomen mee te bouwen. Beide problemen kunnen aangepakt worden door gebruik te maken van een backbone, gebaseerd op de Open Tree of Life. We kunnen dan namelijk de totale data in stukjes knippen die corresponderen met de voornaamste groepen in OpenTOL, om zo per groep een boom te bouwen die we dan op de backbone plakken, en we kunnen de OpenTOL vervolgens ook als een constraint gebruiken voor het algoritme waarmee we bomen bouwen. De opdracht is om dit te prototyperen voor de COI-5P merker in dieren. Conceptueel komt de benadering overeen met een rapport van studenten van LIACS. Een implementatie op basis van data uit BOLD en OpenTOL zou fungeren als pilot
 voor een project-in-voorbereiding voor biodiversa+.
 
## Scripts

### [retrieve_bold_data.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/retrieve_bold_data.py)
Queries the [BOLD Systems API](http://v3.boldsystems.org/index.php/resources/api?type=webservices#combined) to retrieve sequence data for DNA barcodes.  BOLD species names are mapped to [Open Tree Of Life](https://github.com/OpenTreeOfLife/germinator/wiki/Open-Tree-of-Life-Web-APIs#tnrs-methods) taxonomy. FASTA files of the barcodes are made per family group. These are saved as 'fasta/family/{family_name}.fasta'
