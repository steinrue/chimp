
import numpy as np
import sys
import re

np.random.seed(101)

def main():
    '''
    in_file= "test_data3.vcf"
    out_file= "test_data3_pseudo.vcf"

    pseudo(in_file, out_file)

    
    line="1	2226	.	A	T	.	PASS	.	GT	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	1|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|1	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0	0|0"

    oline = ""
    pieces = re.split( '\s',line)

    oline = oline + pieces.pop(0)
    for piece in pieces:
        oline = oline + "\t"
        diploid = re.search('.[|].',piece)
        if(diploid):
            haps = re.split('[|]', diploid.group(0))
            pseudohap = haps[np.random.randint(0,2)]
            oline = oline + pseudohap
            print(str(haps) + " --> " + pseudohap)
        else:
            oline = oline + piece

    print(oline)
    print(line)
    print(line==oline)

    '''
        

def pseudo(in_file, out_file):
    with open(out_file,"w") as ofile:
        with open(in_file, "rt") as ifile:
            line = ifile.readline()
            while(line):

                # parse line for diploids, create pseudohaploid, and write to outfile
                oline = ""
                pieces = re.split( '[\t]',line)

                count_even = True
                head_bool = False
                oline = oline + pieces.pop(0)

                for piece in pieces:
                    if count_even:
                        oline = oline + "\t"
                    else:
                        if head_bool:
                            oline=oline + "-"
                        else:
                            oline = oline + "|"
                        
                    diploid = re.search('.[|].',piece)
                    if(diploid): # if this piece of string is x|y where x,y is {0,1}
                        haps = re.split('[|]', diploid.group(0))
                        pseudohap = haps[np.random.randint(0,2)]
                        oline = oline + pseudohap
                        count_even = not count_even
                    elif(re.search('tsk',piece)):
                        oline = oline + piece
                        head_bool=True
                        count_even = not count_even
                    else:
                        if(piece != ""):
                            oline = oline + piece

                ofile.write(oline)
                line = ifile.readline()


   
if __name__=="__main__":
    
    # name of folder where we want to read files and create consolidated result
    if len(sys.argv) == 1:
        main()
    else:
        in_file = sys.argv[1]
        out_file = sys.argv[2]
        pseudo(in_file, out_file)

        
    





    





    


