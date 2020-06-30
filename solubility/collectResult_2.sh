echo "aqua"
tail $(grep "Overall" ./random_aqua_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
tail $(grep "Overall" ./scaffold_aqua_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2

echo "phys"
tail $(grep "Overall" ./random_phys_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
tail $(grep "Overall" ./scaffold_phys_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2

echo "esol"
tail $(grep "Overall" ./random_esol_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
tail $(grep "Overall" ./scaffold_esol_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2

echo "ochem"
tail $(grep "Overall" ./random_ochem_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
tail $(grep "Overall" ./scaffold_ochem_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2

echo "aqsol"
tail $(grep "Overall" ./random_aqsol_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
tail $(grep "Overall" ./scaffold_aqsol_1/*.newout* | grep "cure" | sort -k 5  | head -n 1 | cut -d":" -f 1) | grep Spearman | cut -d"(" -f 2 | cut -d"," -f 1 | grep -v Spearman | cut -d"=" -f 2
