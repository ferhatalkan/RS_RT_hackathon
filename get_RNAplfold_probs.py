
#!/usr/bin/env python

# The script that parses .ps file of RNAplfold output and reports them in a packed format of your choice
import sys, argparse, re

def RNAplfold_parser(dpfile):
    lines = []
    with open(dpfile) as in_f:
        lines = in_f.readlines()
    i = 0
    seqlen = None
    bps = []
    while i<len(lines):
        line = lines[i]
        if (str(line[1:8])=="winSize") and str(line.rstrip().split()[2])=="def":
            seqlen = int(line.split()[1])
            for k in range(seqlen):
                bps.append([0]*seqlen)

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


def main():
    parser = argparse.ArgumentParser(description='Parse RNAplfold _dp.ps output.')
    parser.add_argument('-r','--RNAplfold_dp_ps_file', type=str, help='RNAplfold result file to parse from.', required=True)
    #parser.add_argument('-b','--pos_specific_bp_prob', type=str, help='Store total base-pair probability of every position to the given output file. (Space seperated values in a single line)')
    #parser.add_argument('-s','--structuredness', type=str, help="Store structuredness of every position to the given output file. (Space seperated values (Total, 5'bps, 3'bps) in each row (for every position))")
    args = parser.parse_args()

    print(RNAplfold_parser(args.RNAplfold_dp_ps_file))

# Run the main when executed
if __name__ == "__main__":
    main()
