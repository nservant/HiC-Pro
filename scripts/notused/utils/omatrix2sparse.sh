
awk '
NR > 1 {
   nr_1 = NR-1;
   for (ii = 2; ii <= NF; ii++) {
      if ($ii && ii-1>=nr_1) {
         printf("%d\t%d\t%d\n", nr_1, ii-1, $ii);
      }
   }
}'

