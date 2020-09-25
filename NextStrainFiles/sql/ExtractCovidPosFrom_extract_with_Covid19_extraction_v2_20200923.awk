#!/usr/bin/awk -f

#Exemple de commande
#./ExtractCovidPosFrom_extract_with_Covid19_extraction_v2_20200923.awk  extract_with_Covid19_extraction_v2_20200923.txt > extract_with_Covid19_extraction_v2_20200923_CovidPos.txt


BEGIN{FS="\t";OFS="\t"}{if(length($19)==0){$19="0"}};{if($15 ~ /^Détecté$|^détecté$/ || $16 ~ /^Détecté$|^détecté$/ || $17 ~ /^Détecté$|^détecté$/ || $18 !~ /NO RESULT/){print $0}}

