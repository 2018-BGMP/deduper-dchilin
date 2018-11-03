########### python script ###############################################################################################

######## arg parse #####################################################################################################

import argparse
def get_arguments():
    parser = argparse.ArgumentParser(description="Demultiplex")
    parser.add_argument('-f', '--input_file', help='insert absolute file path', required=True, type=str)
    parser.add_argument('-o', '--output_file', help='insert absolute file path', required=False, type=str)
    parser.add_argument("-p", "--paired_end", help="designates file is paired end(not single-end)", required=False, type=str)
    parser.add_argument("-s", "--single_end", help="designates file is single end(not paired-end)", required=False, type=str)
    parser.add_argument("-u", "--umi_list", help="designates file containing the list of UMIs", required=False, type=str)
    return parser.parse_args()

#code when using arg parse in bash:
#$ ./script_name -f /path/filename -o other_arguments


############## variables ###################################################################################################
#variables needed when parsing the file
args = get_arguments()
input_file = args.input_file
output_file = args.output_file
umi_list = args.umi_list
p = args.paired_end
s = args.single_end

#creat a variable to call your file 
#since we are using arg parse, you dont need to do this (just comment it out for now)
#f = '/projects/bgmp/dchilin/Bi624/deduper/files/sorted/resorted/sorted_1.sam'

#umi file (its just commented out for now)
#u='/projects/bgmp/dchilin/Bi624/deduper/files/sorted/resorted/UMI.txt'

#create/initiate empty variables to keep track of current & previous chromosome 
current_chr = ''
previous_chr = ''

########### Dictionaries ######################################################################################################

#create dictionary to store all your real UMIs that were given to you 
#u is defined on top via arg parse

UMI_dict={}
for line in open(umi_list).readlines():
    umi=line.strip()
    UMI_dict[umi]=None


#dictionary to check if there are any duplicates
#key: QNAME, Chrom, and POS
dup_dict = {}


################ Functions ###########################################################################################

#some functions require RegEx, so import that pacakge
import re
#isolating UMI from QNAME
#example :  NS500451:154:HWKTMBGXX:1:11101:25533:1187:GTTCACCT
#  UMI : GTTCACCT
def isolate_UMI(string):
    #UMI is found at the very end of the string, use $
    UMI = re.findall(r'([A-Z]+$)', string)[0]
    #re.findall returns all matches as a list. 
    #return [0], even thouhg there is only one things in the list 
    return (UMI)


#checking for single/paired ends
#bitwise flag 1 means strand is paired-end
def check_read(flag):
    read= ((int(flag )& 1) == 1)
    #if True: paired-end (flag equals 1)
    #if False: single-end (flag does not equal 1)
    return read
           
#check if strand is forward or reverse 
#bitwise flag 16 is reverse strand 
def check_strand_read(flag):
    strand = ((int(flag) & 16) == 16)
    #if True: Reverse Strand (flag equals 16)
    #if False: Forward Strand (flag does not equal 16)
    return strand     
    
     
#check cigar string for forward strand ONLY
def check_cigar_forward_string(string):
    #position start at the left most end
    #softclipping can only occur on either ends of the string
    #only account for S thats on the left side of the string
    #ignore any S found on the the right side of the string (it will not affect the position)
    #N, I, & D can only happen in the middle of the read, so dont worry too much about these here 
        #it will afect the position
    #find S's before any M's or any other letter
    S = re.findall(r'(^\d+[S])', string)
    #isolate the number before S
    S = sum([ int(i.strip('S'))  for i in S])
    return S


def check_cigar_reverse_string(string):
    #position start on left most of the string, however this is a reverse string 
    #you want to adjust position to were read is ending, which is the right side
    #find S's on the right side of the string.
    #Ignore the S's on the left (they do not affect the position of the string)
    #find M, D, & N
    S = re.findall(r'(\d+[S]$)', string)
    M = re.findall(r'(\d+[M])', string)
    D = re.findall(r'(\d+[D])', string)
    N = re.findall(r'(\d+[N])', string)
    #isolate the numbers
    S = sum([ int(i.strip('S'))  for i in S])
    M = sum([ int(i.strip('M'))  for i in M])
    D = sum([ int(i.strip('D'))  for i in D])
    N = sum([ int(i.strip('N'))  for i in N])
    #sum S, M, D, N
    total = (S + M + D + N)
    return total
    #note: we do not account for any insertions (I) because it does not consume(affect) the reference read
   

########### line counters #####################################################################################
#line loop is going through 
#we just need this for the first line the script is reading through to set previuos_chr to something)
ctr = 1

#for undetermined/bad UMI's
bad_umi = 0

#for paired_end reads
paired_reads = 0

#for reverse reads
reverse_read = 0

#for forwards reads 
forward_read = 0

#single reads we want to keep
single_reads = 0

#for duplicated reads we dont want to keep
duplicated_reads = 0

################ code ##############################################################################################
#open same file
#[:-4] is used to removed the last 4 characters in the file (.sam)
with open (input_file, 'r') as fh, \
    open(output_file + "_deduped.sam", "w+") as wh, \
    open(output_file + '_unwanted_duped_reads.sam', 'w+') as d_reads, \
    open(output_file + '_undertermined_UMIs.txt', 'w+') as b_umi:
    i = 0
    for line in fh:
        i += 1
        #ignore header lines (start with @)
        #spit those out to deduped file
        if line.startswith('@'):
            wh.write(line)
        else:
            #separate columns from each line into indexes by tab delimination
            #remember index starts with 0
            column = line.split('\t')
            #create variables for each index
            # we only want QNAME, FLAG, RNAME, POS, and CIGAR
            #the others are not needed, but i just labeled it for future reference
            QNAME = column[0] #UMI
            FLAG = column[1] #single/pair, +/-
            RNAME = column[2] #chromosome
            POS = column[3] #position
            MAPQ = column[4]
            CIGAR = column[5] #soft clipping
            RNEXT = column[6]
            PNEXT = column[7]
            TLEN = column[8]
            SEQ = column[9]
            QUAL = column[10]
			            
            ### memory efficiency ##################################################################################
            
            #this is for memory efficiency
            #idea is go through all of one chromosome, once you reach a different chromosome, empty the dictionary
            #set current chromosome to RNAME (index that has the chromosome)
            current_chr = RNAME
            #set previous chromose to the current chromosome
            #this only happens with the very first line reaches this part of the code
            #this is a way to start current and previous chromosome the first chromoomse on the file
            if ctr == 1: 
                previous_chr = current_chr
            ctr += 1
            #once you go through one chromosome, and reach a new chromosome we want to empty the dictioanry
            if current_chr != previous_chr:
                #if current chromosome does not equal previous chromosome(meaning we have encounter a new chromosome)
                #we will set our previous chromosome to the current one
                previous_chr = current_chr
                #emtpy dictionary
                dup_dict = {}

            ### UMI validation ########################################################################################  
            
            #isolate the UMI from QNAME  using the function for it
            #then verify if UMI exists in the UMI dictionary 
            UMI = isolate_UMI(QNAME)
            if UMI not in UMI_dict:
                #UMI not in dictionary, its no good
                #ignore this UMI, throw it to a dump file, and keep track of bad umis
                b_umi.write(line)
                bad_umi += 1
            else:
                 #UMI exists in the dicitonary, continue with code

                ### Paired vs Single-end reads #########################################################################
                
                #check if read is paired-end/single end
                if check_read(FLAG) == True:
                #if true, this is paired-end read
                #code cant handle paired reads, print an error message and quit
                    print('ERROR: Paired-end read')
                    paired_reads += 1
                    break
                else:
                #if not true, this is single-end read, continue with code

                    ###Forward/Reverse strands #########################################################################
                    ### (adjusting postion if needed ) #################################################################

                    #check if strand is Forward/Reverse
					if check_strand_read(FLAG) == True:
                    #if FLAG is 16, its a reverse strand
						reverse_read += 1
						cigar = check_cigar_reverse_string(CIGAR)
                        #adjPOS = int(POS) - 1
						adjPOS = int(POS) + int(cigar) - 1
                        #add '-' to adjPOS to differentiate from the forward strand in the dictionary 
						adjPOS = str(adjPOS) + '-'
					else:
                    #if not 16, its a forward strand
						forward_read += 1
						if 'S' in CIGAR:
							S = check_cigar_forward_string(CIGAR)
                            #subtract the S from POS so that it can adjust for those numbers that were "clipped out"
							adjPOS = int(POS) - int(S)
							#indicate adjPOS that this is a forward(+) strand (to help differentiate from the reverse strand )
							adjPOS = str(adjPOS) + '+'
						else:
							#set POS equal to a variable so we can call on the cigars that dont have S's
							adjPOS = str(POS) + '+'

					### duplicate verfication ###################################################################################

					#all positions should be adjusted
					#add specific columns to dictionary to check for duplicates
					#Only interested if UMI, flag, chromosome and position matches
					#separate it by commas and '' so its easier to read
					#key = (str(UMI))+','+"'"+str(strand)+','+"'"+str(RNAME)+"'"+','+"'"+str(adjPOS)+"'" 
           			#alternatie quicker format than above
					key = '{},{},{}'.format(UMI, RNAME, str(adjPOS))
					#check if there are any duplicate reads in the duplicate dictionary
					if key not in dup_dict:
               			 #not in dictionary, add it to the dictionary
						dup_dict[key] = 1
                		#for memory efficiency, output this line to the deduper file (this is the file with all the wanted silge reads )
						wh.write(line)
               			 #increment counter
						single_reads +=1
					else:
                 		#keys are in the dictionary, there is a mathc. Get rid of this line.
                 		#we only want to keep one copy of this read.
               			 #throw out this read to the unwanted dupliated reads 
						d_reads.write(line)	
                		#increment counter 
						duplicated_reads +=1

print('undertermined umi:', bad_umi) 
print('paired-end reads:', paired_reads)
print('forward reads:', forward_read)                   
print('reverse reads:', reverse_read)
print('unwanted duplicates:', duplicated_reads)
print('wanted deduped reads:', single_reads)


### out put ######################################################################################################################

#***note : 	wanted + unwanted + undetermined umis + header = total lines in original file

#dataset1: 
#lines in original file :1,013,204
#lines in sorted file : 1,013,204
#header(lines with @) : 24
	# undertermined umi: 8,190
	# paired-end reads: 0
	# forward reads: 1,004,945
	# unwanted duplicates: 369,005
	# wanted deduped reads: 635,985

#dataset 2 : 
#lines in orignal file : 1,382,133
#lines in sorted file : 1,382,133
#lines wit#header(lines with @) : 24
	# undertermined umi: 9,319
	# paired-end reads: 0
	# reverse reads: 53
	# unwanted duplicates: 628,602
	# wanted deduped reads: 744,188	

#dataset 3 :
#lines in orignial file: 5,721,174
##lines in sorted file: 5,721,174
#header(lines with @) : 24
	# undertermined umi: 48,939
	# paired-end reads: 0
	# reverse reads: 13,876
	# unwanted duplicates: 1,618,567
	# wanted deduped reads: 4,053,644