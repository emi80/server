# RNA API test dataset

`encodeEvaluation2ga.py` is a simple script to get the RNA API test dataset. It converts a subset of the ENCODE evaluation dataset to a format understandable by the GA4GH reference server.

The script will take care of downloading the data tarball and extract the required information from there. Just run it as follows to get the results in your current folder:

```shell
$ ./encodeEvaluation2ga.py
== Downloading tarball from http://genome.crg.es/~epalumbo/ENCODE-benchmark-data.tgz
== Writing rnaseq tables
== Writing counts tables
== Writing gene expression tables
== DONE
```
