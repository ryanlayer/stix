#!/usr/bin/env python3

import argparse
import os
import gzip
import sys
from collections import Counter, OrderedDict

parser = argparse.ArgumentParser(prog="Convert VCF to query table for STIX")
parser.add_argument("-i", "--input", help="Input VCF file", required=True)
parser.add_argument("-o", "--output", help="Output file", required=True)
parser.add_argument("-c", "--comment", help="Comments appended in ID",default=None,type=str)
parser.add_argument(
    "--pad", help="Padding bases for the left and right breakpoint for all SVTYPE except INS", default=10, type=int)
parser.add_argument(
    "--allow_no_pass", help="Allow variants with out PASS flag", action='store_true')
parser.add_argument("--allow_imprecise",
                    help="Allow IMPRECISE", action='store_true')
parser.add_argument("--allow_emptyGT",
                    help="Allow GTs such as ./.", action='store_true')
parser.add_argument("--allow_all",
                    help="Allow all variants", action='store_true')
# parser.add_argument("--allow_allSVTYPE", help="Allow all SVTYPEs such as complex SVs which are not in <DEL>,<INS>,<DUP>,<INV>,<BND>",action='store_true')
parser.add_argument(
    "--id", help="Which info fields should be used in the ID column, use ',' link multiple", default="")

args = parser.parse_args()

assert os.path.exists(args.input), f"Not found input file:{args.input}"
record_idx = 0
input_records_count = 0
success_records_count = []
inf = None
if os.path.splitext(args.input)[1] == ".gz":
    inf = gzip.open(args.input, 'rt')
else:
    inf = open(args.input)


with open(args.output, 'w') as outf:
    for line in inf:
        record_idx += 1
        xline = line.strip()
        if xline.startswith("#"):
            continue
        xline_sp = xline.split("\t")
        input_records_count +=1
        if len(xline_sp) < 10:
            print(
                "Ignore line at {record_idx} with less than 10 columns",
                file=sys.stderr)
            continue

        chro = xline_sp[0]
        pos = xline_sp[1]
        ID = xline_sp[2]
        refseq = xline_sp[3]
        altseq = xline_sp[4]
        flag = xline_sp[6]
        info_str = xline_sp[7]
        info_list_raw = info_str.split(";")
        info_dict = {}
        for term in info_list_raw:
            if "=" in term:
                k, v = term.split("=")
                info_dict[k] = v
            else:
                info_dict[term] = True
        format_str = xline_sp[8]
        GTidx = 0
        for x in format_str.split(":"):
            if x.strip().upper() != "GT":
                GTidx += 1
                continue
            else:
                break
        sample_str = xline_sp[9]
        sample_list = sample_str.split(":")
        GT = sample_list[GTidx]

        SVTYPE = info_dict.get("SVTYPE", None)
        if SVTYPE is None:
            print(f"Ignore record with missing SVTYPE at line {record_idx}",
                  file=sys.stderr)

        outline = OrderedDict()
        outline['lchr'] = chro
        outline['lstart'] = int(pos) - args.pad
        outline['lend'] = int(pos) + args.pad
        outline['rchr'] = None
        outline['rstart'] = None
        outline['rend'] = None
        outline['svlen'] = None
        outline['svtype'] = SVTYPE

        id_want = [x for x in args.id.strip().split(",") if x.strip() != ""]
        id_out_str = "|".join(
            [f"{k}={info_dict.get(k)}" for k in id_want if info_dict.get(k, "None") is not None])
        outline['ID'] = f"GT={GT}|ID={ID}"
        if len(id_want) > 0:
            outline['ID'] += "|"+ id_out_str
        if args.comment is not None:
            outline['ID'] += "|"+ args.comment

        # for END filed, only SVs which are not BND has this info
        if SVTYPE == "BND":
            clean_second_chrpos = altseq.strip().strip(
                "N").replace("]", "").replace("[", "")
            second_chr = clean_second_chrpos.split(':')[0]
            second_pos = int(clean_second_chrpos.split(':')[1])
            outline['rchr'] = second_chr
            outline['rstart'] = second_pos - args.pad
            outline['rend'] = second_pos + args.pad
            outline['svlen'] = 0
            pass
        elif SVTYPE in ['DEL', 'INV', "DUP"]:
            END = info_dict.get("END", None)
            SVLEN = info_dict.get("SVLEN", None)
            if END is not None:
                outline['rchr'] = chro
                outline['rstart'] = int(END) - args.pad
                outline['rend'] = int(END) + args.pad
                if SVLEN is not None:
                    outline['svlen'] = abs(int(SVLEN))
                else:
                    outline['svlen'] = abs(int(pos) - int(END))
            else:
                if SVLEN is not None:
                    outline['rchr'] = chro
                    outline['rstart'] = int(pos) + abs(int(SVLEN)) - args.pad
                    outline['rend'] =  int(pos) + abs(int(SVLEN)) + args.pad
                    outline['svlen'] = abs(int(SVLEN))
                else:
                    print(f"Ignore record without END and SVLEN at same time at line {record_idx}",
                          file=sys.stderr)

            pass
        elif SVTYPE == "INS":
            SVLEN = info_dict.get("SVLEN", None)
            if SVLEN is not None:
                outline['lchr'] = chro
                outline['lstart'] = int(pos) - args.pad * 2
                outline['lend'] = int(pos)
                outline['rchr'] = chro
                outline['rstart'] = int(pos)
                outline['rend'] = int(pos) + abs(int(SVLEN))
                outline['svlen'] = abs(int(SVLEN))
            else:
                outline['lchr'] = chro
                outline['lstart'] = int(pos) - args.pad * 2
                outline['lend'] = int(pos)
                outline['rchr'] = chro
                outline['rstart'] = int(pos)
                outline['rend'] = int(pos)
                outline['svlen'] = 0
                print(f"Processed insertion without SVLEN at line {record_idx}",
                      file=sys.stderr)

            pass
        # elif args.allow_allSVTYPE:

        #     pass
        else:
            print(f"Ignore record with complex SVTYPE and 'allow_allSVTYPE' set as False at line {record_idx}",
                  file=sys.stderr)
            continue
        
        if args.allow_all != True:

            if args.allow_no_pass != True:
                if flag != "PASS":
                    print(f"Ignore record with non PASS at line {record_idx}",
                        file=sys.stderr)
                    continue

            if args.allow_emptyGT != True:
                if GT not in ['0/0', '1/1', '0/1', '1/0']:
                    print(f"Ignore record with wrong GT:{GT} at line {record_idx}",
                        file=sys.stderr)
                    continue

            if args.allow_imprecise != True:
                if info_dict.get("IMPRECISE", None) is not None:
                    print(f"Ignore record IMPRECISE tag at line {record_idx}",
                        file=sys.stderr)
                    continue

        #stats

        success_records_count.append(SVTYPE)

        # output
        outline_str = None
        if outline['svtype'] in ['DEL', 'INV', "DUP"]:
            if outline['lstart'] < outline['rstart']:
                outf.write(
                    f"{outline['lchr']}:{outline['lstart']}-{outline['lend']}\t{outline['rchr']}:{outline['rstart']}-{outline['rend']}\t{outline['svlen']}\t{outline['svtype']}\t{outline['ID']}\n")

            else:
                outf.write(
                    f"{outline['rchr']}:{outline['rstart']}-{outline['rend']}\t{outline['lchr']}:{outline['lstart']}-{outline['lend']}\t{outline['svlen']}\t{outline['svtype']}\t{outline['ID']}\n")
        else:
            outf.write(
                f"{outline['lchr']}:{outline['lstart']}-{outline['lend']}\t{outline['rchr']}:{outline['rstart']}-{outline['rend']}\t{outline['svlen']}\t{outline['svtype']}\t{outline['ID']}\n")


success_records_count_category = Counter(success_records_count)
print(f"All:{input_records_count},After filteration:{sum(list(success_records_count_category.values()))},By SVTYPE:{success_records_count_category}")