Unpacked DarwinCore archive of the Dutch Species Registry. Downloaded on
2019-01-15 from http://api.biodiversitydata.nl/v2/taxon/dwca/getDataSet/nsr

The key file is Taxa.txt, which was filtered through the following procedure:

- import into sqlite
- filtered to retain only species, with an accepted name, and 1a occurrence status

The query that was run was as follows:

```{sql}
select * from Taxa where taxonomicStatus='accepted name' and taxonRank='species' and occurrenceStatus like '1a %';
```

The resulting file was then exported as a CSV file, which is the Taxa.filterd.txt file in this folder.