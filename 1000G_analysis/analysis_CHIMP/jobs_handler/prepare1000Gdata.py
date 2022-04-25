import numpy
import pandas

numpy.random.seed (4711)

print ("===== CMDs TO DOWNLOAD DATA =====")
# download unphased chromosome 1 from NYGC grch38 data
# unphased is too big, so we download the phased one
# and then we scramble the datasets
print ("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz")
vcfFile = "../../../../../../1000genomes/NYGC_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz"

# also download the corresponding index file
print ("wget http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi")
indexFile = "../../../../../../1000genomes/NYGC_grch38/CCDG_14151_B01_GRM_WGS_2020-08-05_chr1.filtered.shapeit2-duohmm-phased.vcf.gz.tbi"

# download file that maps IDs of individuals to populations from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38
print ("MANUAL: download file that maps IDs of individuals to populations from https://www.internationalgenome.org/data-portal/data-collection/30x-grch38")
idFile = "igsr-1000 genomes 30x on grch38.tsv"

# load id file
idTable = pandas.read_csv (idFile, sep='\t')

# extract ids of given populations
popsToSelect = ["FIN", "LWK", "JPT"]

print ("===== BCFTOOLS CMDs TO EXTRACT POPULATION SPECIFIC DATA =====")
for pop in popsToSelect:
    ids = list((idTable.loc[idTable["Population code"] == pop])["Sample name"])
    vcfOutWithInfo = f"info_phased_chr1_{pop}.vcf"
    bcfToolsCmd = f'bcftools view --force-samples -s {",".join(ids)} {vcfFile} > {vcfOutWithInfo}'
    # bcfToolsCmd = f'bcftools view --force-samples -s {",".join(ids)} -r chr2:{int(130e6)}-{int(140e6)} {vcfFile} > phased_chr1_{pop}.vcf'
    print (bcfToolsCmd)
    vcfOutWithNoSeg = f"noseg_phased_chr1_{pop}.vcf"
    removeInfoCmd = f'bcftools annotate -x INFO {vcfOutWithInfo} > {vcfOutWithNoSeg}'
    print (removeInfoCmd)
    vcfOutOnlySeg = f"phased_chr1_{pop}.vcf"
    removeNonSegCmd = f'egrep "0[/|][1-9]|[1-9][/|]0|[1-9][/|][1-9]|#" {vcfOutWithNoSeg} > {vcfOutOnlySeg}'
    print (removeNonSegCmd)
    vcfOutUnphased = f"unphased_chr1_{pop}.vcf"
    unphaseCmd = f"python unphase.py {vcfOutOnlySeg} {vcfOutUnphased} {numpy.random.randint (999999999)}"
    print (unphaseCmd)
    print (pop, len(ids))

