#!/usr/bin/env python
import gzip
import FileUtilities as fu

ACCEPTED_CHR = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20","21","22", "X", "Y", "MT"]

def count_alt(strng):
    strng = strng.upper()
    lst = list(strng)
    a=0
    c=0
    t=0
    g=0
    n=0
    ast=0
    match_sum=0
    vc='SNP'

    for l in lst:
        l=str(l)
        if l=='A':
            a=a+1

        elif l=='C':
            c=c+1

        elif l=='T':
            t=t+1

        elif l=='G':
            g=g+1

        elif l=='*':
            ast=ast+1
            
        elif l=='N':
            n=n+1

        elif str(l)=='.' or str(l)==',':
            match_sum=match_sum+1

        elif l=='+' or l=='-':
            vc='INDEL'

    return str(match_sum) + '\t' + str(a) + '\t' + str(c) + '\t' + str(t) + '\t' + str(g) + '\t' + str(ast) + '\t' + str(n) + '\t' + str(vc)


def get_count_alt_qpileup(pileup):
    outfile=pileup+'.alt_count'
    fh = gzip.open(pileup+'.gz', 'r')
    fu.delete(outfile+'.gz')
    fh_out = gzip.open(outfile+'.gz', 'w')
    header = 'CHROM' + '\t' + 'POS' + '\t' + 'REF' + '\t' + 'DEPTH' + '\t' + 'REFCOUNT' + '\t' + 'A' + '\t' + 'C' + '\t' + 'T' + '\t' + 'G' + '\t' + 'DEL' + '\t' + 'N' + '\t' + 'VC'
    fh_out.write(header+'\n')

    for line in fh:
        fields=line.split('\t')
        #1	909768	A	41	GGgGgggGGGGGGggggGGggggGggggggggggggggg^]g^]g
        tofile=str(fields[0])+'\t'+str(fields[1])+'\t'+str(fields[2])+'\t'+str(fields[3])+'\t'+count_alt(str(fields[4]))
        fh_out.write(tofile+'\n')

    fh.close()
    fh_out.close()



####Find the maximum allele(s) out of A,T,C,G,DEL
def find_max_alleles(allele_counts, allele_list):
    a=allele_counts
    b=allele_list
    max_a=max(a)
    max_b=""
    for i in range(5):
        if (a[i]==max_a):
            max_b=max_b+b[i]
            
    return max_b


def get_macount(altcount):
    outfile_aap = altcount + '.aap'
    outfile_ref = altcount + '.ref'
    fh = gzip.open(altcount+'.gz', 'r')

    fu.delete(outfile_aap+'.gz')
    fu.delete(outfile_ref+'.gz')

    fh_aapout = gzip.open(outfile_aap+'.gz', 'w')
    fh_refout = gzip.open(outfile_ref + '.gz', 'w')

    # in "alt_count file"
    #CHROM POS REF DEPTH REF COUNT A C T G DEL N VC
    header='CHROM'+'\t'+'POS'+'\t'+'REF'+'\t'+'DEPTH'+ '\t' + 'REFCOUNT' + '\t' + 'A' + '\t' + 'C' + '\t' + 'T' + '\t' + 'G' + '\t' + 'DEL' + '\t' + 'N' + '\t' + 'VC' + '\t' + 'AD' + '\t' + 'MACOUNT' + '\t' + 'MA'
    fh_aapout.write(header+'\n')
    fh_refout.write(header+'\n')

    linecount = 0
    alt_alleles=['A', 'C', 'T', 'G', 'DEL']
    
    for line in fh:
        line = str(line).strip()
        if linecount > 0:
            fields=line.split('\t')
            ichr=str(fields[0])
            if ichr.startswith('chr'):
                chr=ichr[3:len(ichr)]
            else:
                chr=ichr

            pos=str(fields[1])
            ref=str(fields[2])
            depth=int(fields[3])
            ref_count=int(fields[4])
            n_count=int(fields[10])
            vc=str(fields[11])
            ## reference is not counted in the previous step. In pileup, if bp==ref then it is '.' or ','
            alleles_counts = map(int, fields[5:10]) #fields[5:10]=>A,C,T,G,DEL
            sum_of_alternative=sum(alleles_counts)
                
            if ( fu.find_first_index(ACCEPTED_CHR, chr.strip())>-1 ):
                # Determine "max_allele", "max_allele_count", and "adjusted_depth"
                if (sum_of_alternative == 0):
                    adjusted_depth = int(depth - n_count)
                    max_allele_count=0
                    max_allele="."
                    tofile=line+'\t'+str(adjusted_depth)+'\t'+str(max_allele_count)+'\t'+str(max_allele)
                    fh_refout.write(tofile+'\n')
                            
                elif (sum_of_alternative > 0):
                    max_allele_count = max(alleles_counts) # maximum alternative alele count
                    #max_allele = alt_alleles[alleles_counts.index(max_allele_count)] # maximum allele
                    max_allele=find_max_alleles(alleles_counts, alt_alleles)
                    adjusted_depth = int(ref_count + max_allele_count) #adjusted depth
                    tofile=line+'\t'+str(adjusted_depth)+'\t'+str(max_allele_count)+'\t'+str(max_allele)
                    fh_aapout.write(tofile+'\n')

        linecount = linecount+1
    
    fh.close()
    fh_aapout.close()
    fh_refout.close()

