import argparse

def get_arguments():
    parser = argparse.ArgumentParser(description="Python script to remove PCR duplicates from SAM file.")
    parser.add_argument("-f", help="set filename", required=True, type=str)
    parser.add_argument("-p", help="set single or paired end reads", required=False, type=bool)
    parser.add_argument("-u", help="set UMIs", required=False, type=str)
    #parser.add_argument("-h", help="type for help", required=False, type=str)

    return parser.parse_args()







def bit_checker(bit):
    '''Takes the bitwise flag and checks it for strandedness. Assumes read is mapped,
    and data are single-end. Returns "+" or "-" depending on strand.'''

    #Check if read is mapped
    if (bit & 4) == 4:
        raise ValueError("Read is unmapped")

    strand = "+"

    #Changes strand value if bit 16 is used
    if (bit & 16) == 16:
        strand = "-"

    return strand






def cigar_parse(cigar):
    '''Takes CIGAR string and returns the amount of soft clipping present. Only checks for left side soft clipping
    and soft clipping less then 10.'''
    soft_clip = 0

    #Check if soft clipping present and set amount
    if cigar[1] == 'S':
        soft_clip = cigar[0]
    return int(soft_clip)



def get_UMI(qname):
    '''Parse out UMI from the QNAME field of the SAM file'''
    qname = qname.strip()
    fields = qname.split(":")
    return fields[7]

def erase_dict(u):
    '''Reset values of the keys to empty'''
    dict = {}
    with open(u)as uf:
        for umi in uf:
            umi = umi.strip()
            dict[umi] = []
    return dict



args = get_arguments()

#initialize variables and data structures
filename = args.f
paired = args.p
u = args.u

if paired == True:
    raise Exception("Paired end files not yet implemented.")

UMI_dict = erase_dict(u)

current_chr = ""
bad_umi = 0
dup = 0

#open SAM file
with open(filename)as fh:
    #open output file
    with open(filename + "_deduped", "w")as output:
        #split SAM file into lines

        for line in fh:
            #skip header lines
            if line[0] != '@':
                #remove newline chars
                line = line.strip()

                #split lines into columns
                parts = line.split("\t")

                #get data from columns
                chr = parts[2]

                #if new chr reached reset dictionary keys to have empty values
                if chr != current_chr:
                    current_chr = chr
                    UMI_dict = erase_dict(u)
                    

                #get corrected pos and strand
                position = str(int(parts[3]) - cigar_parse(parts[5]))
                strand = str(bit_checker(int(parts[1])))
                pair = [position, strand]




                #get UMI from QNAME
                UMI = get_UMI(parts[0])


                #check for unknown UMIs
                if UMI in UMI_dict.keys():
                    #check if duplicate
                    if pair in UMI_dict[UMI]:
                        # count duplicate
                        dup = dup + 1
                    else:
                        # if not duplicate add pair to dict and output line
                        UMI_dict[UMI].append(pair)
                        output.write(line + "\n")
                else:
                    #count unknown UMI
                    bad_umi = bad_umi + 1


print("duplicates: " + str(dup))
print("bad UMIs: " + str(bad_umi))