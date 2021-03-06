#!/usr/bin/env bash

#expands a bedgraph file to single nucleotide intervals of file $2
# and uses bedtools map to count up features over file $1
# using sum over column $3
# output goes to file $4
#ie line from file $2:
#   chrI 10 12 20 2 +
#will be converted to:
#   chrI 10 11 20 2 +
#   chrI 11 12 20 2 +
# and counts from this file over interval in $1
#   chrI 5 20 . . +
# will give chrI 5 20 . . + 4

awk '{OFS="\t"}
{ start=$2;
  end=$3;
  for(i=start;i<end;i++){
    $2=i;
    $3=i+1;
    print $0
  }
}' $2 | \
sort -k1,1 -k2,2n | \
bedtools map -a $1 -b - -c $3 -o sum > $4
