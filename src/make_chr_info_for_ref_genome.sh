#!/bin/bash

echo -e "chr\tlength\n\
1\t249250621\n\
2\t243199373\n\
3\t198022430\n\
4\t191154276\n\
5\t180915260\n\
6\t171115067\n\
7\t159138663\n\
8\t146364022\n\
9\t141213431\n\
10\t135534747\n\
11\t135006516\n\
12\t133851895\n\
13\t115169878\n\
14\t107349540\n\
15\t102531392\n\
16\t90354753\n\
17\t81195210\n\
18\t78077248\n\
19\t59128983\n\
20\t63025520\n\
21\t48129895\n\
22\t51304566" > chr_length_tmp.txt

{
echo -e "chr\tlength\tcumulative_sum"

awk -F"\t" 'BEGIN {sum=0} NR>1 {print $0 "\t" sum; sum += $2}' chr_length_tmp.txt
} > chr_length.txt
