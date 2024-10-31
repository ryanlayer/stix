from collections import OrderedDict
import argparse
from copy import deepcopy
import os

VERSION="1.0.0"

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i', "--input", nargs="*",
                        help="input individual STIX annotated vcfs")

    parser.add_argument('-o', "--output",
                        help="output file")
    parser.add_argument('-v', '--version', 
                    action='version',
                    version=f'%(prog)s {VERSION}')
    
    args = parser.parse_args()

    if len(args.input) < 2:
        print("please provide at least two vcf files to merge!")
        exit(1)

    if os.path.exists(args.output):
        answer = input(f"output file:{args.output} exists! override?[Yy/Nn]")
        if answer.strip().lower() == "y":
            os.remove(args.output)
        else:
            print("Nothing to do,exit...")
            exit(0)

    for single_vcf in args.input:
        if not single_vcf.strip().upper().endswith("VCF"):
            print(f"Warning: detected wrong suffix in {single_vcf}")

    merged_vcf_content = OrderedDict()
    merged_vcf_header = []
    clean_input = [x for x in args.input if x != args.output]
    base_vcf = clean_input[0]
    

    # process the first input sharded vcf file as template vcf.
    count = 0
    with open(base_vcf) as baseinf:
        for line in baseinf:
            count += 1
            if line.startswith('#'):
                merged_vcf_header.append(line)
                continue

            xline = line.strip().split("\t")
            uniKey = "@".join(xline[0:5])
            info_str = xline[7].split(";")
            info_dict = {x.split("=")[0]: x.split("=")[1]
                         for x in info_str if '=' in x}

            STIX_ONE = info_dict.get("STIX_ONE", None)
            STIX_ZERO = info_dict.get("STIX_ZERO", None)
            if (STIX_ONE is not None) and (STIX_ZERO is not None):
                merged_vcf_content[uniKey] = {
                    "raw_xline": deepcopy(xline),
                    "stix_count": [int(STIX_ZERO), int(STIX_ONE)]
                }
            else:
                raise Exception(
                    f"Can not get STIX count from base vcf:{base_vcf},at line {count}")

    # iter read different vcfs  and extract the STIX_ZERO and STIX_ONEs
    count = 0
    for other_vcf in clean_input[1:]:
        with open(other_vcf) as otherinf:
            for line in otherinf:
                count += 1
                if line.startswith("#"):
                    continue
                xline = line.strip().split("\t")
                uniKey = "@".join(xline[0:5])
                info_str = xline[7].split(";")
                info_dict = {x.split("=")[0]: x.split("=")[1]
                             for x in info_str if '=' in x}
                STIX_ONE = info_dict.get("STIX_ONE", None)
                STIX_ZERO = info_dict.get("STIX_ZERO", None)
                # print(STIX_ZERO,STIX_ONE)
                if (STIX_ONE is not None) and (STIX_ZERO is not None):
                    # print("xx")
                    merged_vcf_content[uniKey]["stix_count"][0] += int(STIX_ZERO)
                    merged_vcf_content[uniKey]["stix_count"][1] += int(STIX_ONE)
                else:
                    raise Exception(
                        f"Can not get STIX count from base vcf:{base_vcf},at line {count}")


    # iter read different vcfs and extract the STIX_ZERO and STIX_ONEs
    with open(args.output, 'w') as outf:
        no_stix_freq = True
        for ee in merged_vcf_header:
            if "STIX_FREQ" in ee:
                no_stix_freq = False
        if no_stix_freq:
            merged_vcf_header.insert(-1,'##INFO=<ID=STIX_FREQ,Number=1,Type=Float,Description="Population frequency from STIX annotation. equals to STIX_ONE/(STIX_ONE+STIX_ZERO)">\n')
        for x in merged_vcf_header:
            outf.write(x)
            outf.flush()
        for y in merged_vcf_content.values():
            raw_list = y['raw_xline']
            stix_count = y["stix_count"]
            raw_info_list = raw_list[7].split(";")
            # print(raw_info_list)
            # print(raw_info_list,stix_count)
            new_info_list = []

            for each_info in raw_info_list:
                if "STIX_ONE" in each_info:
                    new_info_list.append(f"STIX_ONE={stix_count[1]}")
                    continue
                elif "STIX_ZERO" in each_info:
                    new_info_list.append(f"STIX_ZERO={stix_count[0]}")
                    continue
                elif "STIX_FREQ" in each_info:
                    continue # do nothing
                else:
                    new_info_list.append(each_info)
            new_info_list.append(f"STIX_FREQ={stix_count[1]/(stix_count[1]+stix_count[0]) :.4f}") 
            raw_list[7] = ";".join(new_info_list)
            out_str = '\t'.join(raw_list) + "\n"
            outf.write(out_str)


if __name__ == "__main__":
    main()
