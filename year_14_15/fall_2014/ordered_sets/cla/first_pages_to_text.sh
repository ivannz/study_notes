#! /usr/bin/env bash
cd ./CLA/cla.inf.upol.cz/papers
## Find all PDFs and execute pdfttotext (from xpdf) to extract
##  the first page of the document and save into a uniquely
##  named but easily idetifiable txt file
find ./*/*.pdf -type f -exec pdftotext -f 1 -l 1 -enc UTF-8 -eol unix -nopgbrk {} {}.txt \;
## Store all generated TXTs into a TAR archive so that the
##   directory structure is preserved.
tar -cf ./abstracts.tar ./*/*.txt
## Remove the now useless first pages
rm -rf ./*/*.txt
## Copy the TAR into a more convenient location (the original working directory)
cp ./abstracts.tar ./../../../
cd ./../../../
## Extract the TAR contents into a special directory
mkdir ./abstracts
cd ./abstracts
tar -xf ./../abstracts.tar
cd ./..
## Find all lines starting from "Abstract" (case-sensitive) and store
##  the found lines into a file
ggrep -no -m 1 -P 'Abstract.*' ./abstracts/*/*.txt > ./abstracts.txt
## Cut the found lines into a tabulated file
perl -i -pe 's/^\.\/abstracts(.*?\.pdf)\.txt:(\d+):\s*[Aa]bstract[\.:]?\s+(.*?)$/cla.inf.upol.cz\/papers\1\t\2\t\3/g' ./abstracts.txt
## Do the smae for keywords
ggrep -no -m 1 -P 'Keywords.*' ./abstracts/*/*.txt > ./Keywords.txt
perl -i -pe 's/^\.\/abstracts(.*?\.pdf)\.txt:(\d+):\s*[Kk]eywords[\.:]?\s+(.*?)$/cla.inf.upol.cz\/papers\1\t\2\t\3/g' ./Keywords.txt
#### grep -noi abstract ./abstracts/*/*.txt | wc -l
#### ^(.*?\.pdf)\.txt:(\d+):\s*[Aa]bstract[\.:]?\s+(.*?)$
## List all the papers' PDFs
ls ./*/*.pdf > ./../../../all_cla_pdf.txt
