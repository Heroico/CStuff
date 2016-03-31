# 1000 Genomes download

First, download data via FTP
```bash
# Connect from a sane environment with access to internet such as syncmon
# FTP download
$ ftp ftp://ftp.1000genomes.ebi.ac.uk
ftp> cd vol1/ftp/release/20130502/
ftp> mget *
```

Then you can convert by:

```bash
$ plink --vcf chr1.vcf --out chr1
```