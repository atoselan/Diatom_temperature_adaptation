module load bcftools


# Fixed het positions. I.e. 0/1 to 0/0 or 1/1
while read sample; do 
  bcftools view \
    --types snps evoExp.bcfCall.vars.bcf | \
  bcftools filter -e 'GT=="./."' | \
  bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS [%GT ] [%AD{0} %AD{1} ]\n' | \
  awk '{
    d1 = $7 + $8
    d2 = $9 + $10
    d3 = $11 + $12
    d4 = $13 + $14
    if (d1 >= 10 && d2 >= 10 && d3 >= 10 && d4 >= 10) {
      if ((($3=="0/1" && ($8/d1>=0.1)) &&
           ($4=="0/1" && ($10/d2>=0.1)) &&
           ($5=="0/1" && ($12/d3>=0.1))) &&
          (($6=="0/0" && ($13/d4>=0.99)) ||
           ($6=="1/1" && ($14/d4>=0.99))))
        print $0
    }
  }' | wc -l
done < TP_one_sample.txt


# 0/0 in T0 to 0/1 or 1/1 in evolved cell line
while read sample; do
  bcftools view \
    --types snps evoExp.bcfCall.vars.bcf | \
  bcftools filter -e 'GT=="./."' | \
  bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS [%GT ] [%AD{0} %AD{1} ]\n' | \
    awk '{
      d1 = $7 + $8
      d2 = $9 + $10
      d3 = $11 + $12
      d4 = $13 + $14
      if (d1 >= 10 && d2 >= 10 && d3 >= 10 && d4 >= 10) {
        if ((($3=="0/0" && ($7/d1>=0.99)) &&
             ($4=="0/0" && ($9/d2>=0.99)) &&
             ($5=="0/0" && ($11/d3>=0.99))) &&
            (($6=="0/1" && ($14/d4>=0.1)) ||
             ($6=="1/1" && ($14/d4>=0.99))))
          print $0
      }
    }' | wc -l
done < TP_sample_names.txt


# 1/1 to 0/0 or 0/1 in evolved cell line
while read sample; do
  bcftools view \
    --types snps evoExp.bcfCall.vars.bcf | \
  bcftools filter -e 'GT=="./."' | \
  bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS [%GT ] [%AD{0} %AD{1} ]\n' | \
    awk '{
      d1 = $7 + $8
      d2 = $9 + $10
      d3 = $11 + $12
      d4 = $13 + $14
      if (d1 >= 10 && d2 >= 10 && d3 >= 10 && d4 >= 10) {
        if ((($3=="1/1" && ($8/d1 >= 0.99)) &&
             ($4=="1/1" && ($10/d2 >= 0.99)) &&
             ($5=="1/1" && ($12/d3 >= 0.99))) &&
            (($6=="0/1" && ($14/d4 >= 0.1)) ||
             ($6=="0/0" && ($13/d4 >= 0.99))))
          print $0
      }
    }'  | wc -l
done < TP_sample_names.txt


# Check for tri-allelic snps 0/0 to 0/2 or 1/2 or 2/2
while read sample; do
  bcftools view \
    --types snps \
    evoExp.bcfCall.vars.bcf | \
  bcftools filter -e 'GT=="./."' | \
  bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS [%GT ] [%AD{0} %AD{1} %AD{2}]\n' | \
  awk 'NF >= 18 {
    d1 = $7 + $8 + $9
    d2 = $10 + $11 + $12
    d3 = $13 + $14 + $15
    d4 = $16 + $17 + $18
    if (d1 >= 10 && d2 >= 10 && d3 >= 10 && d4 >= 10) {
      if (($3=="0/0" && $7/d1 >= 0.99) &&
          ($4=="0/0" && $10/d2 >= 0.99) &&
          ($5=="0/0" && $13/d3 >= 0.99)) {
        if (($6=="0/2" && $18/d4 >= 0.1) || ($6=="1/2" && $18/d4 >= 0.1) || ($6=="2/2" && $18/d4 >= 0.99)) {
          print $0
        }
      }
    }
  }' | wc -l
done < TP_sample_names.txt


# Novel allelic variants
while read sample; do
  bcftools view \
    --types snps \
    evoExp.bcfCall.vars.bcf | \
  bcftools filter -e 'GT=="./."' | \
  bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS [%GT ] [%AD{0} %AD{1} %AD{2}]\n' | \
  awk 'NF >= 18 {
    d1 = $7 + $8 + $9
    d2 = $10 + $11 + $12
    d3 = $13 + $14 + $15
    d4 = $16 + $17 + $18
    if (d1 >= 10 && d2 >= 10 && d3 >= 10 && d4 >= 10) {
      # T0 samples must all be 0/0 with allele 0 >=99%
      if (($3 == "0/0" && $7 / d1 >= 0.99) &&
          ($4 == "0/0" && $10 / d2 >= 0.99) &&
          ($5 == "0/0" && $13 / d3 >= 0.99)) {
        
        # Evolved sample must be 0/2, 1/2, or 2/2 with allele 2 ≥10%
        if ($6 == "0/2" && $18 / d4 >= 0.1 || 
            $6 == "1/2" && $18 / d4 >= 0.1 || 
            $6 == "2/2" && $18 / d4 >= 0.99) {
          
          # Check for novel SNP gain (0/1 or 1/1 to 0/2, 1/2, or 2/2)
          if (($3 == "0/1" || $3 == "1/1") && ($6 == "0/2" || $6 == "1/2" || $6 == "2/2")) {
            print $0
          }
        }
      }
    }
  }' | wc -l
done < TP_sample_names.txt


# Get deletions
while read sample
do
  bcftools view -i 'TYPE="indel" && strlen(ALT) < strlen(REF)' evoExp.bcfCall.vars.bcf | \
    bcftools filter -e 'GT=="./."' | \
    bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS %REF %ALT [%GT ] [%AD{0} %AD{1} ]\n' | \
    sed 's/\./0/g' | \
    awk '{ if (($8 != $5 && $8 != $6 && $8 != $7) && (($15+$16>=10) && $16/($15+$16)>0.1))  print $0 }' | wc -l
done < TP_sample_names.txt


# Get insertions
while read sample
do
  bcftools view -i 'TYPE="indel" && strlen(ALT) > strlen(REF)' evoExp.bcfCall.vars.bcf | \
    bcftools filter -e 'GT=="./."' | \
    bcftools query \
    --samples T0_22C_1,T0_22C_2,T0_22C_3,${sample} \
    --format '%CHROM %POS %REF %ALT [%GT ] [%AD{0} %AD{1} ]\n' | \
    sed 's/\./0/g' | \
    awk '{ if (($8 != $5 && $8 != $6 && $8 != $7) && (($15+$16>=10) && $16/($15+$16)>0.1))  print $0 }' | wc -l
done < TP_sample_names.txt
