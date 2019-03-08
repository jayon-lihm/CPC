import os
import sys
import FileUtilities as fu
import functions_alt_count as AltCnt

### Tabulate pileup files (list and count all allternative alleles from pileup files)
## INPUT: chr${ichar}_q${MapQual}_d${depth}_mapped.pileup
## OUTPUT: chr${ichar}_q${MapQual}_d${depth}_mapped.pileup.alt_count
q_d_pileup=sys.argv[1]
AltCnt.get_count_alt_qpileup(pileup=q_d_pileup)

