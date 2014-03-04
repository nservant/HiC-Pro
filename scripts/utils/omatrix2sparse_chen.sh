
awk '
{
   for (ii = 1; ii <= NF; ii++) {
      if ($ii) {
         printf("%d\t%d\t%d\n", NR, ii, $ii);
      }
   }
}'

