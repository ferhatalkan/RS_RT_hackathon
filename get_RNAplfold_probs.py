
#!/usr/bin/env python

# The script that parses .ps file of RNAplfold output and reports them in a packed format of your choice
import sys, argparse, re, os
from statistics import mean

def RNAplfold_parser(dpfile):
    lines = []
    with open(dpfile) as in_f:
        lines = in_f.readlines()
    i = 0
    seq = ''
    bps = []
    while i<len(lines):
        line = lines[i]
        if line[1:9]=="sequence":
            while lines[i+1].startswith(") } def")==False:
                i += 1
                line = lines[i]
                seq += re.sub('\\\n','',line.replace('\\','')).translate(str.maketrans('ACGTUacgtu','ACGTTACGTT'))

            for k in range(len(seq)):
                bps.append([0]*len(seq))

        elif lines[i].startswith("%start of base pair probability data"):
            while lines[i+1].startswith("showpage")==False:
                i += 1
                cols = lines[i].split()
                if len(cols)<4:
                    sys.stderr.write("#ERROR: RNAplfold base pair row does not have enough(3) columns.")
                    print(cols)
                    exit()
                else:
                    bps[int(cols[0])-1][int(cols[1])-1] = float(cols[2])
                    bps[int(cols[1])-1][int(cols[0])-1] = float(cols[2])
        i += 1
    return bps

def get_str_feature(BPmatrix, CDSs, CDSe, w=50, maxW=200, inside=False):
    result = []
    seqlen = len(BPmatrix)

    allsums = []
    for i in range(0, min((CDSe+100),seqlen)-1):
        if inside == False:
            allsums.append(sum(BPmatrix[i][i:min(i+maxW, seqlen)]))
        elif inside == True:
            allsums.append(sum(BPmatrix[i][i:min(i+w, i+maxW, seqlen)]))

    for cod_i in range(int((CDSe-CDSs)/3)):
        #print(cod_i,seqlen, len(allsums), CDSs, CDSe, min(CDSs+(cod_i*3)+16, seqlen-1),  min(CDSs+(cod_i*3)+15+w, seqlen))
        result.append(str(mean([allsums[nt_i] for nt_i in range(min(CDSs+(cod_i*3)+16, seqlen-2), min(CDSs+(cod_i*3)+15+w, seqlen-1))])))
    return(result)

def main():
    parser = argparse.ArgumentParser(description='Parse RNAplfold _dp.ps output.')
    parser.add_argument('-r','--RNAplfold_dp_ps_folder', type=str, help='RNAplfold result folder to parse from.', required=True)
    #parser.add_argument('-b','--pos_specific_bp_prob', type=str, help='Store total base-pair probability of every position to the given output file. (Space seperated values in a single line)')
    #parser.add_argument('-s','--structuredness', type=str, help="Store structuredness of every position to the given output file. (Space seperated values (Total, 5'bps, 3'bps) in each row (for every position))")
    parser.add_argument('-o','--output_f', type=str, help='Output file path.', required=True)
    args = parser.parse_args()

    with open(args.output_f,'w') as outf:
        all_files = os.listdir(args.RNAplfold_dp_ps_folder)
        for fi in range(len(all_files)):
            if fi % 100 == 0:
                print(fi, "out of", len(all_files))
            RNAplfold_file = all_files[fi]
            cols = RNAplfold_file.split("_")
            if cols[4].endswith('001'):
                #print(RNAplfold_file)
                SM = RNAplfold_parser(args.RNAplfold_dp_ps_folder+"/"+RNAplfold_file)
                for i in range(len(cols)):
                    if cols[i] == "CDS":
                        CDSs, CDSe = [int(n) for n in cols[i+1].split("-")]
                        break
                outf.write("\t".join([cols[0]]+get_str_feature(SM,CDSs,CDSe))+'\n')

# Run the main when executed
if __name__ == "__main__":
    main()
