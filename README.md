# barcode-constrained-phylogeny
This repository contains code and data for building topologically-constrained phylogenies.

NEEDS TRANSLATING! Also need to add the two files to repo?

De internationale database BOLD heeft DNA barcodes voor honderdduizenden soorten (met vaak meerdere barcodes per soort). 
Als je de data in BOLD zou filteren op merker, er zitten immers meerdere merkers in de database, en die vervolgens per merker zou alignen, dan zou dit redelijke data zijn om fylogenetische bomen mee te bouwen. Echter, er zijn twee beperkingen: 1) er zit een bovengrens aan de grootte van alignments en bomen die te bouwen zijn, en 2) de barcodes zijn vrij kort dus geven vrij weinig signaal om grote bomen mee te bouwen. Beide problemen kunnen aangepakt worden door gebruik te maken van een backbone, gebaseerd op de Open Tree of Life. We kunnen dan namelijk de totale data in stukjes knippen die corresponderen met de voornaamste groepen in OpenTOL, om zo per groep een boom te bouwen die we dan op de backbone plakken, en we kunnen de OpenTOL vervolgens ook als een constraint gebruiken voor het algoritme waarmee we bomen bouwen. De opdracht is om dit te prototyperen voor de COI-5P merker in dieren. Conceptueel komt de benadering overeen met een rapport van studenten van LIACS. Een implementatie op basis van data uit BOLD en OpenTOL zou fungeren als pilot
 voor een project-in-voorbereiding voor biodiversa+.
 
## Scripts
### [unzip_targz.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/unzip_targz.py)
Unzips a targz file, more specifically a [datarelease](https://www.boldsystems.org/index.php/datapackage?id=BOLD_Public.30-Dec-2022) from BOLD Systems containing a snapshot of the barcode database (more than 8 million barcodes as of 30-DEC-2022).

### [bold_data_dump.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/bold_data_dump.py)
Puts relevent BOLD data columns into a custom SQLite database.


### [alter_tables.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/alter_tables.py)
Manipulates the BOLD data in the custom database and makes two tables, one for taxon data and one for barcode entries.


### [map_to_opentol.py](https://github.com/naturalis/barcode-constrained-phylogeny/blob/main/src/map_to_opentol.py)
Uses the [Checklisbank API](https://api.checklistbank.org/) to map BOLD taxon names to [Open Tree of Life](https://tree.opentreeoflife.org/opentree/argus/opentree13.4@ott93302) taxonomy IDs. 
