#!/usr/bin/env bash

#expands a bedgraph file to single nucleotide intervals of file $2
# and uses bedtools map to count up features over file $1
# using column $3
#ie line:
#   chrI 10 12 20 x +
#will be converted to:
#   chrI 10 11 20 x +
#   chrI 11 12 20 x +

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
