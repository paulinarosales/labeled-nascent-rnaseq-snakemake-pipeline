from pysam import FastaFile
import argparse
import re

parser = argparse.ArgumentParser(description="RAC filter for m6A data")

parser.add_argument('--bed_file', type=str, help="Bed file to be filtered")
parser.add_argument('--fasta_file', type=str, help="Fasta file with genome")
parser.add_argument('--window_size', type=int, default= 0, help="How big the window should be. Use 0 for regular strict RAC filter.")
parser.add_argument('--adapt_coordinates', type=str, default=True, help="Whether we should switch coordinates to A in RAC motif")
parser.add_argument('--single_cell', type=str, default=False, help="If we want to have an additional barcode column")



"""
Slightly simplified version of the windowed RAC filter
Does not add coverage and dart mutations, just filters for RAC in a window and adapts coordinates
"""




args = parser.parse_args()

window_size = args.window_size
print_header = True
# Print out rac motif A coordinates instead
adapt_coordinates = args.adapt_coordinates
single_cell = args.single_cell


def correct_coordinates(seq,strand,testing=False):
    # Find the RAC motif closest to transition
    if strand=="-":
        possible_motifs = [m.start() for m in re.finditer("GTT", seq)]+[m.start() for m in re.finditer("GTC", seq)]
    else:
        possible_motifs = [m.start() for m in re.finditer("AAC", seq)]+[m.start() for m in re.finditer("GAC", seq)]
    # Determine position of closest RAC motif to deaminated C/G
    pos_in_seq = min(possible_motifs,key=lambda x:abs(x-window_size))
    correction_factor = pos_in_seq-window_size
    # Scale to get the A in the RAC motif instead of the motif beginning
    correction_factor +=1
    if testing:
        print(possible_motifs)
        print(pos_in_seq)
        print(correction_factor)
        print(seq[window_size+correction_factor])
    if strand == "+":
        assert seq[window_size+correction_factor]=="A", "Error. No adenosine found!"
    else:
        assert seq[window_size+correction_factor]=="T", "Error. No adenosine found!"
    return correction_factor



fasta = FastaFile(args.fasta_file)

# Bed file needs at least 6 columns

if print_header:
    # We only keep certain columns we are interested in
    if single_cell:
        print("#chr"+"\t"+"start"+"\t"+"end"+"\t"+"gene_region_etc"+"\t"+"edit_ratio"+"\t"+"strand"+"\t"+"old_site_id"+"\t"+"barcode")
    else:
        print("#chr"+"\t"+"start"+"\t"+"end"+"\t"+"gene_region_etc"+"\t"+"edit_ratio"+"\t"+"strand"+"\t"+"old_site_id")
    
    
    #print("chr"+"\t"+"start"+"\t"+"end"+"\t"+"gene_region_etc"+"\t"+"edit_ratio"+"\t"+"strand"+"\t"+"control_edit_ratio"
    #      +"\t"+"control_coverage"+"\t"+"dart_edit_ratio"+"\t"+"dart_coverage"+"\t"+"experiment"+"\t"+"control_mutations"+"\t"+"dart_mutations")

# Site IDs from meyer can be non-unique, so we make our own ones
    
    
    
    
# Iterate over bed file line by line
with open(args.bed_file,"r") as f:
    for line in f.readlines():
        line = line.strip()

        # Case when header line is still in the file
        if line[:4]=="#chr":
            continue
            
        # IMPORTANT: When window size is 0, this means, we will expand only ot the 3nts
        # This means, for - strand we have to reverse complement
        if window_size == 0:
            if line.split("\t")[5]=="+":
                start = int(line.split("\t")[1])-2
                end = int(line.split("\t")[2])
            else:
                start = int(line.split("\t")[1])
                end = int(line.split("\t")[2])+2
                
            seq = fasta.fetch(line.split("\t")[0],start,end).upper()
        else:    
            start = int(line.split("\t")[1])-window_size
            # Case where chr end is already reached
            if start < 0:
                start = 0
            # Get fasta sequence (expand 5nts in each direction
            seq = fasta.fetch(line.split("\t")[0],start,int(line.split("\t")[2])+window_size).upper()
        # Find RAC motif in genomic sequence
        if line.split("\t")[5]=="+":
            if "AAC" in seq or "GAC" in seq:
                site_id = str(line.split("\t")[0])+"_"+str(line.split("\t")[1])+"_"+str(line.split("\t")[2]+"_"+str(line.split("\t")[5]))
                #dart_mutations = round(float(line.split("\t")[8])*float(line.split("\t")[9]))
                #control_mutations = round(float(line.split("\t")[6])*float(line.split("\t")[7]))
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                # Printing the methylated A instead of the editing
                if adapt_coordinates==True:
                    correction_factor = correct_coordinates(seq,line.split("\t")[5])
                    start = str(int(start)-correction_factor)
                    end = str(int(end)-correction_factor)
                #else:
                #    start = str(int(line.split("\t")[1]))
                #    end = str(int(line.split("\t")[2]))
                    
                if single_cell:
                    barcode = line.split("\t")[11]
                    print(line.split("\t")[0]+"\t"+start+"\t"+end+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]
                          +"\t"+line.split("\t")[5]+"\t"+site_id+"\t"+barcode)
                else:
                    print(line.split("\t")[0]+"\t"+start+"\t"+end+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]
                          +"\t"+line.split("\t")[5]+"\t"+site_id)

           
        

        # Look for reverse (complement!) if on minus strand
        # C-Ts on reverse strand are processed as G-As (see IGV tracks), so we have to look for reverse-complemented motif
        else:
            if "GTC" in seq or "GTT" in seq:
                site_id = str(line.split("\t")[0])+"_"+str(line.split("\t")[1])+"_"+str(line.split("\t")[2]+"_"+str(line.split("\t")[5]))
                #dart_mutations = round(float(line.split("\t")[8])*float(line.split("\t")[9]))
                #control_mutations = round(float(line.split("\t")[6])*float(line.split("\t")[7]))
                start = line.split("\t")[1]
                end = line.split("\t")[2]
                # Printing the methylated A instead of the editing
                if adapt_coordinates==True:
                    correction_factor = correct_coordinates(seq,line.split("\t")[5])
                    start = str(int(start)+correction_factor)
                    end = str(int(end)+correction_factor)
                if single_cell:
                    barcode = line.split("\t")[11]
                    print(line.split("\t")[0]+"\t"+start+"\t"+end+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]
                          +"\t"+line.split("\t")[5]+"\t"+site_id+"\t"+barcode)
                else:
                    print(line.split("\t")[0]+"\t"+start+"\t"+end+"\t"+line.split("\t")[3]+"\t"+line.split("\t")[4]
                          +"\t"+line.split("\t")[5]+"\t"+site_id)
           
       
