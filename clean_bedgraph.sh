#!/bin/sh

#averages overlapping intervals that may exist in a bedgraph file
#useful for example to make a lifted bedgraph file amenable for bedGraphToBigWig ucsc tool
#input bedgraph = $1; outputfile=$2

awk '{OFS="\t"}{
  if((c==$1 && $2>=e) || c!=$1){
    if(use_ar == 1){
      si=s;
      ei=s+1;
      vi=ar_v[s];
      ni=ar_n[s];
      #print "    doing interval"c"\t"s"\t"e;
      #for(i=s;i<e;i++){
      #  print c"\t"i"\t"i+1"\t"ar_v[i]" ./. "ar_n[i];
      #}
      for(i=s;i<e;i++){
        if(ar_v[i]==vi && ar_n[i]==ni){
          ei=i+1;
        }else{
          print c"\t"si"\t"ei"\t"vi/ni;
          si=ei;ei=si+1;vi=ar_v[i];ni=ar_n[i];
        }
      }
      print c"\t"si"\t"ei"\t"vi/ni;
    }else{
      print c"\t"s"\t"e"\t"v;
    }
    c=$1;s=$2;e=$3;v=$4;
    #print "   x1"c"\t"s"\t"e"\t"v;
    delete ar_v;
    delete ar_n;
    use_ar = 0;
  }else{
    #print "use_ar: "use_ar
    if(use_ar == 0){
      #print "   adding first range"c"\t"s"\t"e
      for(i=s;i<e;i++){
        ar_v[i] = v;
        ar_n[i] = 1;
        #print "   adding first"c"\t"i"\t"ar_v[i]"\t"ar_n[i];
      }
      use_ar = 1;
    }
    for(i=$2;i<$3;i++){
      ar_v[i] += $4;
      ar_n[i] += 1;
      #print "   appending first"c"\t"i"\t"ar_v[i]"\t"ar_n[i];
    }
    e=$3;
  }
}END{
  if(use_ar == 1){
      si=s;
      ei=s+1;
      vi=ar_v[s];
      ni=ar_n[s];
      for(i=s;i<e;i++){
        if(ar_v[i]==vi && ar_n[i]==ni){
          ei=i+1;
        }else{
          print c"\t"si"\t"ei"\t"vi/ni;
          si=ei;ei=si+1;vi=ar_v[i];ni=ar_n[i];
        }
      }
      print c"\t"si"\t"ei"\t"vi/ni;
    }else{
      print c"\t"s"\t"e"\t"v;
    }
}' $1 | sed 1d > $2


