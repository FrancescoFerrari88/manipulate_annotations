#!/Users/ferrari/anaconda3/bin/python
# see that happens
import argparse
import sys
import os
import textwrap


def parse_args(defaults={"verbose":False,
                         "feature_to_extract":"genes",
                         "from_where":"gene",
                         "before":0,
                         "after":0,
                         "out_file":None,
                         "percentile_range":None}):

    """parse arguments from the command line"""

    parser = argparse.ArgumentParser(
        prog=sys.argv[0],
        formatter_class=argparse.RawDescriptionHelpFormatter,
        #description=textwrap.dedent(__description__),
        add_help= False
    )


    # positional - required
    parser.add_argument("gtf_file",
                       metavar="GTF",
                       help="gencode gtf file from which you want to extract genomic coordinates")
    # general arguments
    general = parser.add_argument_group('general arguments')
    general.add_argument("-o", "--output_file_name",
                         dest="out_file",
                         help="output file name",
                         default=defaults["out_file"])
    general.add_argument("-h","--help",
                         action="help",
                         help="show this help message and exit")
    general.add_argument("-v","--verbose",
                         dest="verbose",
                         action="store_true",
                         help="verbose output (default: {})".format(defaults["verbose"]),
                         default=defaults["verbose"])
    general.add_argument("-f","--feature",
                         dest="FEATURE",
                         choices=["genes","TSS","TES"],
                         default=defaults["feature_to_extract"],
                         help="the feature that you want to extract form the gtf file")
    general.add_argument("-w","--from_what",
                         dest="from_what",
                         choices=["gene","transcript"],
                         default=defaults["from_where"],
                         help="Do you want to extract TSS/TES info from genes or transcripts?")
    general.add_argument("-b","--before_feature",
                         dest="BEFORE_FEATURE",
                         type=int,
                         help="number of bp to include before the feature of interest",
                         default=defaults["before"]
                         )
    general.add_argument("-a","--after_feature",
                         dest="AFTER_FEATURE",
                         type=int,
                         help="number of bp to include after the feature of interest",
                         default=defaults["after"]
                         )
    return parser





def info_to_dict(info):

    dict_info={}

    lista = info.split(";")
    lista = [i.strip() for i in lista]

    for j in lista:
        pair = j.split()
        if len(pair) == 2:
            key,value = pair[0],pair[1]

            if not key in dict_info:
                dict_info[key] = value
            else:
                pass
                #print("Warning: duplicate keys have been found in the information section: {}".format(key))
        else:
            if len(pair) != 0:
                print("Warning: more than two values found: {}".format(pair))

    return dict_info





def line_to_dict(gtf_line_list):

    keys = ["chromosome_name",
            "annotation_source",
            "feature_type",
            "genomic_start_location",
            "genomic_end_location",
            "score",
            "genomic_strand",
            "genomic_phase_CDS",
            "additional_information"]

    line_dict = {}

    if len(gtf_line_list) == 9:
        for i in range(len(keys)):
            line_dict[keys[i]] = gtf_line_list[i]
        line_dict["additional_information"] = info_to_dict(line_dict["additional_information"])
    else:
        if not gtf_line_list[0].startswith("#"):
            print("Warning! The following line does not conform to the expected gencode gtf format:\n{}".format(" ".join(gtf_line_list)))

    return line_dict




def make_bed(list_dict, feature, from_where, before, after):

    bed_out_list = []

    if feature == "genes" and list_dict["feature_type"] == "gene":
        if list_dict["genomic_strand"] == '+':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - before),
                str(int(list_dict["genomic_end_location"]) + after),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - after),
                str(int(list_dict["genomic_end_location"]) + before),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

    elif feature == "TSS" and list_dict["feature_type"] == from_where:
        if list_dict["genomic_strand"] == '+':
             bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - before),
                str(int(list_dict["genomic_start_location"]) + after),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_end_location"]) - after),
                str(int(list_dict["genomic_end_location"]) + before),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

    elif feature == "TES" and list_dict["feature_type"] == from_where:
        if list_dict["genomic_strand"] == '+':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_end_location"]) - before),
                str(int(list_dict["genomic_end_location"]) + after),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

        elif list_dict["genomic_strand"] == '-':
            bed_out_list = [
                list_dict["chromosome_name"],
                str(int(list_dict["genomic_start_location"]) - after),
                str(int(list_dict["genomic_start_location"]) + before),
                list_dict["additional_information"]["gene_id"],
                list_dict["score"],
                list_dict["genomic_strand"]
            ]

    return bed_out_list

def main():

    parser = parse_args()
    arg = parser.parse_args()
    print(arg.gtf_file, arg.FEATURE)

    if arg.out_file == None:
        gtf_file_name = "_".join(os.path.basename(arg.gtf_file).split(".")[:-1])

        arg.out_file = "{}_{}_{}_window{}-{}.bed".format(gtf_file_name, arg.FEATURE, arg.from_what,arg.BEFORE_FEATURE, arg.AFTER_FEATURE)
    print(arg.out_file)

    with open(arg.gtf_file) as in_file, open(arg.out_file,"w") as out_file:

        for line in in_file:

            lista = line.strip().split("\t")
            list_dict = line_to_dict(lista)
            if len(list_dict) > 0:
                printing_line = make_bed(list_dict, arg.FEATURE, arg.from_what, arg.BEFORE_FEATURE, arg.AFTER_FEATURE)
                if len(printing_line) > 0:
                    printing_line[-1] = printing_line[-1]+"\n"
                    out_file.write("\t".join(printing_line))


if __name__ == "__main__" :
    main()
