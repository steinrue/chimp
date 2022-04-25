import sys
import numpy
import allel
# print(allel.__version__)
import time
import cyvcf2




def unphase (inVcf, outVcf):

    # read the vcf with scikit-allel, just to get number of snps
    print ("[GET_NR_SNPS]")
    print (f"Reading: {inVcf}")
    startTime = time.perf_counter()
    callset = allel.read_vcf (inVcf)
    print (f"Took {(time.perf_counter() - startTime):.2f} seconds.")

    # no tri-allelic?
    assert (sum(callset["variants/ALT"][:,2] != '') == 0)
    assert (sum(callset["variants/ALT"][:,1] != '') == 0)
    assert (sum(callset["variants/ALT"][:,0] == '') == 0)

    snpsInFile = callset["calldata/GT"].shape[0]
    print (snpsInFile)

    print ("[DONE]")

    print ("[UNPHASE]")
    print (f"File to unphase: {inVcf}")
    print (f"Unphased output written to: {outVcf}")

    # go through the vcf
    vcfIFS = cyvcf2.VCF (inVcf)

    # get some randomness
    numIndividuals = len(vcfIFS.samples)
    randomness = numpy.random.randint (2, size=(numIndividuals, snpsInFile))

    # create a new vcf Writer using the input vcf as a template.
    vcfOFS = cyvcf2.Writer (outVcf, vcfIFS)

    count = 0
    allIdxs = numpy.arange (numIndividuals)

    for v in vcfIFS:

        # see what goes
        # what are the indices to be flipped?
        toFlip = allIdxs[randomness[:,count] == 1]
        for idx in toFlip:
            # flip it
            v.genotypes[idx][0], v.genotypes[idx][1] = v.genotypes[idx][1], v.genotypes[idx][0]

        # make sure we have new genotypes
        v.genotypes = v.genotypes
        # and write it
        vcfOFS.write_record(v)

        # increase count
        count += 1
        if (count % 100000 == 0):
            print (count)

    vcfOFS.close()
    vcfIFS.close()

    print ("[DONE]")


def main():

    # parse parameters
    if (len(sys.argv) != 4):
        print ("usage: python <script_name> <inVcf> <outVcf> <seed>")
    else:
        # the file to unphase
        inputVcf = sys.argv[1]
        # output to be written here
        outputVcf = sys.argv[2]
        # and the seed
        numpy.random.seed (int(sys.argv[3]))

        # and do the unphasing
        unphase (inputVcf, outputVcf)


if __name__ == "__main__":
    main()
