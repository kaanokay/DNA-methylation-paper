#!/usr/bin/env python
# coding: utf-8

import pandas as pd
import numpy
import gzip
import argparse
import sys

# # Mapping table:
# ## genome1_coord, Anchor+expanded_cigar_offset, genome2_coord
def processCigar(cigar):
    '''Expand cigar string
    def processCigar(cigar:<str>):
    Parameters:
        cigar:<str> - Samtool Compatible Cigar String. e.g. 2=
    Output:
        out:<str> - expanded cigar string. e.g.:==
    '''
    out = ''
    N = 0
    for symbol in cigar:
        if symbol in '0123456789':
            N = 10*N + int(symbol)
        else:
            #if (symbol != 'D'):
            if (N == 0):
                out += symbol
            else:
                out += N*symbol
            N = 0
    return out

def combineCigar(cigar):
    '''Compress expanded cigar string
    def combineCigar(cigar:<str>):
    Parameters:
        cigar:<str> - Expanded cigar strings. e.g. ==
    Output:
        out:<str> - Samtool Compatible Cigar String. e.g.:2=
    '''
    cigar = cigar +'$'
    out = ''
    N = 0
    start = 0
    for i in range(1,len(cigar)):
        if cigar[i-1] != cigar[i]:
            out += str(i-start) + cigar[i-1]
            start = i
    return out 

def offsetexchange(cigar, ref_i, alt_i, int_i):
    '''Given a compressed cigar string between two sequences, return the base pair mapping between the two sequences
    def offsetexchange(cigar <str>, ref_i <int>, alt_i <int>):
    Parameters:
        cigar: <str> - Samtool Compatible Cigar String. e.g.:2=
        ref_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome 1 (as ref)
        alt_i: <int> - integer, linear coordinates of the first base of the proximal anchor on Genome 2 (as alt)
    Output:
        offsets:<list> - a list of identical coordinates, 
        (1) linear coord on ref, 
        (2) offset on the intermediate string, 
        (3) linear coord on alt)
    '''
    offsets = []
    cigar = processCigar(cigar)
    for i, s in enumerate(cigar): 
        if s == '=':
            offsets.append((ref_i, int_i+i, alt_i))
            alt_i += 1
            ref_i += 1
        
        if s == 'I':
            offsets.append((-1, int_i+i, alt_i))
            alt_i += 1
        if s == 'D':
            offsets.append((ref_i, int_i+i, -1))
            ref_i += 1
        if s == 'X':
            offsets.append((ref_i, int_i+i, alt_i))
            alt_i += 1
            ref_i += 1
            
    return offsets

def load_alignment(alignmentfile):
    '''Load pairwise alignment results between the two genomes
    
    def load_alignment(alignmentfile):
    Parameters:
        alignmentfile <str>: pairwise alignment filename
    Output:
        anchor_dict <dict>: 
        keys: anchor name
        values: cigar strings
    '''
    with open(alignmentfile, 'r') as f:
        data = f.read()
        data = data.split('\n')
        data.pop(0)
    alignment = {}
    for item in data:
        itemlist = item.split(', ')
        if len(itemlist) == 2:
            alignment[itemlist[0]] = itemlist[1]
    return alignment

def get_anchor_intermediate_pos(anchor_pos_g1, anchor_pos_g2, alignment):
    anchorinfo = []
    pos = 1
    for i, (node, pos_g1) in enumerate(anchor_pos_g1[:-1]):
        anchorinfo.append((node, pos))
        a = alignment.get(node, "?")
        if a == "?":
            d_g1 = anchor_pos_g1[i+1, 1] - anchor_pos_g1[i, 1]
            d_g2 = anchor_pos_g2[i+1, 1] - anchor_pos_g2[i, 1]
            d = max(d_g1, d_g2)
            pos += d
            continue
        cigar = processCigar("45=" + a)
        pos += len(cigar)
    
    node = anchor_pos_g1[i+1, 0]
    assert node == "SINK", "Incomplete"
    anchorinfo.append((node, pos))
    anchor_dict = dict(anchorinfo)
    return anchor_dict

def load_anchor(file1, file2, chromo):
    '''Load anchor information for the two genomes
    
    def load_anchor(file1, file2, chromo):
    Parameters:
        file1 <str>: anchor_info filename for g1
        file2 <str>: anchor_info filename for g2
        chromo <str>: current chromosome name
    Output:
        anchor_pos_g1 <array>: anchor name, linear coordinates on Genome A
        anchor_pos_g2 <array>: anchor name, linear coordinates on Genome B
    '''
    anchor_info_g1 = pd.read_csv(file1, index_col = None)
    anchor_info_g2 = pd.read_csv(file2, index_col = None)
    
    anchor_info_g1 = anchor_info_g1[anchor_info_g1['Chromosome'] == chromo]
    anchor_info_g2 = anchor_info_g2[anchor_info_g2['Chromosome'] == chromo]
    
    anchor_pos_g1 = anchor_info_g1[['Anchor', 'Pos']].values
    anchor_pos_g2 = anchor_info_g2[['Anchor', 'Pos']].values
    
    return (anchor_pos_g1, anchor_pos_g2)

def load_position_list(file_pos, chromo, cpg):
    '''Load coordinate list used for mapping
    
    def load_position_list(args):
    Parameters:
        file_pos <str>: input bedMethyl filename
        chromo <str>: current chromosome
    Output:
        poslist <list>: linear coordinates on this chromosome        
    '''
    bedMethyl = pd.read_csv(file_pos, sep = '\t', index_col = None, header = None)
    chromlist = numpy.unique(bedMethyl.iloc[:,0].tolist())
    bedMethyl = bedMethyl[bedMethyl.iloc[:,0] == chromo]
    
    if cpg:
      poslist = bedMethyl.iloc[:,1].tolist()
      return bedMethyl, poslist, chromlist
    else:
      poslist_start = bedMethyl.iloc[:,1].tolist()
      poslist_end = bedMethyl.iloc[:,2].tolist()
      return bedMethyl, poslist_start, poslist_end, chromlist

def construct_mapping_table(alignment, anchor_pos_g1, anchor_pos_g2):
    '''Construct mapping table for the chromosome
    
    def construct_mapping_table(alignment <dict>, anchor_pos_g1 <dict>, anchor_pos_g2 <dict>):
    Parameters:
        alignment <dict>: keys - anchor name; values - cigar string
        anchor_pos_g1 <dict>: keys - anchor name; values - linear coordinates on genome A
        anchor_pos_g2 <dict>:keys - anchor name; values - linear coordinates on genome B
    Output:
        mapping_table <array>: mapping between two genome and the intermediate string
        (col_1) linear coord on ref, 
        (col_2) offset on the intermediate string, 
        (col_3) linear coord on alt        
    '''
    mapping_table = []
    anchorlist = list(alignment.keys())
    anchor_intermediate_dict = get_anchor_intermediate_pos(anchor_pos_g1, anchor_pos_g2, alignment)
    anchor_pos_g1 = dict(anchor_pos_g1)
    anchor_pos_g2 = dict(anchor_pos_g2)
    
    for node in anchorlist:
        ref_coord = anchor_pos_g1[node]
        alt_coord = anchor_pos_g2[node]
        intermediate_coord = anchor_intermediate_dict[node]
        
        if alignment.get(node) == "?":
            continue
        
        cigar = "45=" + alignment[node]
        offsets = offsetexchange(cigar, ref_coord, alt_coord, intermediate_coord)
        mapping_table += offsets
    
    mapping_table = numpy.array(mapping_table)
    mapping_table = mapping_table.astype(int)
    return mapping_table

def mapping_cpg(poslist, mapping_table, col_index, intermediate, reverse):
    '''Given a list of input list, return the output mapping of the coordinates
    
    def mapping_cpg(poslist <list>, mapping_table <array>, anchor_pos <dict>, col_index <int>):
    Parameters:
        poslist <list>: a list of indices that used for mapping
        mapping_table <array>: mapping between two genome and the intermediate string
        (col_1) linear coord on ref, 
        (col_2) offset on the intermediate string, 
        (col_3) linear coord on alt   
        anchor_pos <array>: anchor name; linear coordinates on that genome mapped from
        intermediate<bool>: if True, output intermediate genome coordinates, if false, output the alternative linear genome coordinates
        
    Output:
        list of coordinates on the alternative genome or intermediate genome
    '''
    mapping_poslist = []
    if intermediate:
        if not reverse:
            if col_index == 0:
                mapping_dict = dict(zip(mapping_table[:,0], mapping_table[:,1]))
            else:
                mapping_dict = dict(zip(mapping_table[:,2], mapping_table[:,1]))
        if reverse:
            if col_index == 0:
                mapping_dict = dict(zip(mapping_table[:,1], mapping_table[:,0]))
            else:
                mapping_dict = dict(zip(mapping_table[:,1], mapping_table[:,2]))
    else:
        if col_index == 0:
            mapping_dict = dict(zip(mapping_table[:,0], mapping_table[:,2]))
        else:
            mapping_dict = dict(zip(mapping_table[:,2], mapping_table[:,0]))
    
    for pos in poslist:
        altpos = mapping_dict.get(pos, "")
        mapping_poslist.append((pos, altpos))
        
    return mapping_poslist

def mapping(poslist_start, poslist_end, mapping_table, col_index, intermediate, reverse):
    '''Given a list of input list, return the output mapping of the coordinates
    
    def mapping(poslist <list>, mapping_table <array>, anchor_pos <dict>, col_index <int>):
    Parameters:
        poslist <list>: a list of indices that used for mapping
        mapping_table <array>: mapping between two genome and the intermediate string
        (col_1) linear coord on ref, 
        (col_2) offset on the intermediate string, 
        (col_3) linear coord on alt   
        anchor_pos <array>: anchor name; linear coordinates on that genome mapped from
        intermediate<bool>: if True, output intermediate genome coordinates, if false, output the alternative linear genome coordinates
        
    Output:
        list of coordinates on the alternative genome or intermediate genome
    '''
    mapping_poslist_start = []
    mapping_poslist_end = []
    
    if intermediate:
        if not reverse:
            if col_index == 0:
                mapping_dict = dict(zip(mapping_table[:,0], mapping_table[:,1]))
            else:
                mapping_dict = dict(zip(mapping_table[:,2], mapping_table[:,1]))
        if reverse:
            if col_index == 0:
                mapping_dict = dict(zip(mapping_table[:,1], mapping_table[:,0]))
            else:
                mapping_dict = dict(zip(mapping_table[:,1], mapping_table[:,2]))
    else:  
        if col_index == 0:
            mapping_dict = dict(zip(mapping_table[:,0], mapping_table[:,2]))
        else:
            mapping_dict = dict(zip(mapping_table[:,2], mapping_table[:,0]))
    
    for pos_start in poslist_start:
        altpos_start = mapping_dict.get(pos_start, "")
        mapping_poslist_start.append((pos_start, altpos_start))
        
    for pos_end in poslist_end:
        altpos_end = mapping_dict.get(pos_end, "")
        mapping_poslist_end.append((pos_end, altpos_end))
    
    return mapping_poslist_start, mapping_poslist_end

def update_bedMeth_cpg(bedMeth_unmapped, mapped_poslist):
  assert len(bedMeth_unmapped.index) == len(mapped_poslist)
  
  mapped_poslist_start = [pos[1] for pos in mapped_poslist]
  mapped_poslist_end = [pos[1]+1 if (str(pos[1]) != "" and int(pos[1]) != -1) else "" for pos in mapped_poslist]
  unmapped_poslist = [pos[0] for pos in mapped_poslist if str(pos[1]) == "" or int(pos[1]) == -1]
  
  bedMeth_unmapped.iloc[:,1] = mapped_poslist_start
  bedMeth_unmapped.iloc[:,2] = mapped_poslist_end
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,1] != ""]
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,1] != -1]
  
  assert len(mapped_poslist) == len(bedMeth_unmapped.index) + len(unmapped_poslist)
  
  return bedMeth_unmapped, unmapped_poslist

def update_bedMeth(bedMeth_unmapped, mapped_poslist_start, mapped_poslist_end):
  assert len(bedMeth_unmapped.index) == len(mapped_poslist_start)
  assert len(bedMeth_unmapped.index) == len(mapped_poslist_end)
  
  mapped_poslist_start_out = [pos[1] for pos in mapped_poslist_start]
  mapped_poslist_end_out = [pos[1] for pos in mapped_poslist_end]
  
  unmapped_poslist_start = [pos[0] for pos in mapped_poslist_start if str(pos[1]) == "" or int(pos[1]) == -1]
  unmapped_poslist_end = [pos[0] for pos in mapped_poslist_end if str(pos[1]) == "" or int(pos[1]) == -1]
  unmapped_poslist = unmapped_poslist_start + unmapped_poslist_end
  unmapped_poslist.sort()
  
  bedMeth_unmapped.iloc[:,1] = mapped_poslist_start_out
  bedMeth_unmapped.iloc[:,2] = mapped_poslist_end_out
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,1] != ""]
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,1] != -1]
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,2] != ""]
  bedMeth_unmapped = bedMeth_unmapped[bedMeth_unmapped.iloc[:,2] != -1]
  
  return bedMeth_unmapped, unmapped_poslist

def main():
    # Create parser object
    parser = argparse.ArgumentParser()
 
    # Define arguments for parser object
    parser.add_argument("-b", "--bedMethyl", type = str, nargs = 1,
                        metavar = "file_name", default = None,
                        help = "bedMethyl input file - first 3 columns must be (chr, start, end)", required=True)
     
    parser.add_argument("-p", "--PairwiseAlignment", type = str, nargs = 1,
                        metavar = "path", default = None,
                        help = "Pairwise alignment file", required=True)
    
    parser.add_argument("-t", "--Type", type = str, nargs = 1,
                        metavar = "input coordinate", default = [0],
                        help = "Input coordinate type, 'g1' for genome 1, and 'g2' for genome 2. When Reverse is True, set this to the desired output coordinate type.", required=True)
    
    parser.add_argument("-c", "--Chromosome", type = str, nargs = 1,
                        metavar = "chromosome", default = None,
                        help = "Specify chromosome", required=True)
     
    parser.add_argument("-a", "--Anchor", type = str, nargs = 2,
                        metavar = ('file1','file2'), help = "anchor_info file for genome 1 and genome 2", required=True)
    
    parser.add_argument("-o", "--outputDir", type = str, nargs = 1,
                    metavar = "output directory", help = "Output directory", default = None)
    
    parser.add_argument("-i", "--Intermediate", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), nargs = 1, default=[False],
                    metavar = "output intermediate coordinates", help = "Output intermediate genome coordinates or alternative coordinates, True for outputting intermediate genome coordinates")

    parser.add_argument("-r", "--Reverse", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), nargs = 1, default=[False],
                    metavar = "convert from intermediate coordinates into a reference genome", help = "True for converting from the intermediate coordinates into one of the contributing reference genomes (specified by Type). If Intermediate is False, this does nothing.")
 
    parser.add_argument("-cpg", "--CpG", type=lambda x: (str(x).lower() in ['true', '1', 'yes']), nargs = 1, default=[False],
                    metavar = "mapping CpGs", help = "Are CpGs being mapped? True if mapping CpGs, False if mapping another feature.")

    parser.add_argument("-gn", "--genomeNames", type = str, nargs = 2,
                        metavar = ("g1_name", "g2_name"), default = ['g1', 'g2'], help = "Genome names for g1 and g2")

    args = parser.parse_args()
    
    bedMeth = args.bedMethyl[0]
    pairAlign = args.PairwiseAlignment[0]
    t = args.Type[0]
    chromo = args.Chromosome[0]
    anchor_info_g1 = args.Anchor[0]
    anchor_info_g2 = args.Anchor[1]
    if args.outputDir is None:
      outDir = ""
    else:
      outDir = args.outputDir[0]
    intermediate = args.Intermediate[0]
    reverse = args.Reverse[0]
    g1_name = args.genomeNames[0]
    g2_name = args.genomeNames[1]
    cpg = args.CpG[0]
    
    sys.stdout.write("bedMethyl input file: " + bedMeth + '\n')
    sys.stdout.write("Pairwise Alignment file: " + pairAlign + '\n')
    sys.stdout.write("Anchor info files: " + anchor_info_g1 + ", " + anchor_info_g2 + '\n')
    if not intermediate or not reverse:
      sys.stdout.write("Input bedMethyl reference coordinate system: " + t + '\n')
    if intermediate and reverse:
      sys.stdout.write("Output bedMethyl reference coordinate system: " + t + '\n')
    if (len(outDir) > 0):
      sys.stdout.write("Output directory: " + outDir + '\n')
    else:
      sys.stdout.write("Output directory: Working directory" + '\n')
    sys.stdout.write("Intermediate mapping?: " + str(intermediate) + '\n')
    if intermediate and not reverse:
      sys.stdout.write("Mapping to intermediate?: " + str(intermediate) + '\n')
    if intermediate and reverse:
      sys.stdout.write("Mapping from intermediate?: " + str(intermediate) + '\n')
    sys.stdout.write("Mapping CpGs?: " + str(cpg) + '\n')
    sys.stdout.write("Current chromosome: " + chromo + '\n')
    sys.stdout.write("g1 name: " + g1_name + '\n')
    sys.stdout.write("g2 name: " + g2_name + '\n')
    
    if (len(outDir) > 0):
      if outDir[len(outDir)-1] != "/":
        outDir = outDir + "/"
    
    out_baseName = '.'.join(bedMeth.split("/")[len(bedMeth.split("/"))-1].split(".")[:-1])
    if intermediate:
      if not reverse:
        if t == 'g1':
          outfile = outDir + out_baseName + "_mappedToIntermediateFrom" + g1_name + ".bed"
        elif t == 'g2':
          outfile = outDir + out_baseName + "_mappedToIntermediateFrom" + g2_name + ".bed"
      elif reverse:
        if t == 'g1':
          outfile = outDir + out_baseName + "_mappedFromIntermediateTo" + g1_name + ".bed"
        elif t == 'g2':
          outfile = outDir + out_baseName + "_mappedFromIntermediateTo" + g2_name + ".bed"
    elif t == 'g1' and not intermediate:
      outfile = outDir + out_baseName + "_mappedTo" + g2_name + ".bed"
    elif t == 'g2' and not intermediate:
      outfile = outDir + out_baseName + "_mappedTo" + g1_name + ".bed"
    
    if t == 'g1':
      t = 0
    elif t == 'g2':
      t = 2
    else:
      raise Exception("Invalid input coordinate type. '-t' parameter must be set as 'g1' or 'g2'.")

    anchor_dict = load_alignment(pairAlign)
    
    anchor_pos_g1, anchor_pos_g2 = load_anchor(anchor_info_g1, anchor_info_g2, chromo)
    
    if cpg:
      bedMeth_unmapped, poslist, chromlist = load_position_list(bedMeth, chromo, cpg)
    else:
      bedMeth_unmapped, poslist_start, poslist_end, chromlist = load_position_list(bedMeth, chromo, cpg)
    
    chromindex = chromlist.tolist().index(chromo)
    
    mapping_table = construct_mapping_table(anchor_dict, anchor_pos_g1, anchor_pos_g2)
    
    if cpg:
      mapping_coordinates = mapping_cpg(poslist, mapping_table, t, intermediate, reverse)
      bedMeth_mapped, poslist_unmapped = update_bedMeth_cpg(bedMeth_unmapped, mapping_coordinates)
    else:
      mapping_coordinates_start, mapping_coordinates_end = mapping(poslist_start, poslist_end, mapping_table, t, intermediate, reverse)
      bedMeth_mapped, poslist_unmapped = update_bedMeth(bedMeth_unmapped, mapping_coordinates_start, mapping_coordinates_end)
    
    if chromindex == 0:
      bedMeth_mapped.to_csv(outfile, sep = '\t', mode = 'w', header = False, index = False)
      with open(out_baseName + "_unmappedPositions.tsv", 'w') as fp:
        for pos in poslist_unmapped:
            fp.write(chromo + '\t' + str(pos) + '\n')
    else:
      bedMeth_mapped.to_csv(outfile, sep = '\t', mode = 'a', header = False, index = False)
      with open(out_baseName + "_unmappedPositions.tsv", 'a') as fp:
        for pos in poslist_unmapped:
            fp.write(chromo + '\t' + str(pos) + '\n')
      

if __name__ == "__main__":
    # Calling the main function
    main()
