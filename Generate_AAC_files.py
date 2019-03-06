import os
import sys
import FileUtilities as fu
import functions_alt_count as AltCnt


## Find the largest alternative allele, its counts, adjusted depth (count of ref allele + count of the largest alternative allele)
## INPUT: chr${ichar}_q${MapQual}_d${depth}_mapped.pileup.alt_count  
## OUTPUT: chr${ichar}_q${MapQual}_d${depth}_mapped.pileup.alt_count.aap.gz, chr${ichar}_q${MapQual}_d${depth}_mapped.pileup.alt_count.ref.gz
file_in=sys.argv[1]
AltCnt.get_macount(altcount=file_in)

