#!/bin/bash
rm variant_summary.txt.gz
wget --tries = 100 https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
