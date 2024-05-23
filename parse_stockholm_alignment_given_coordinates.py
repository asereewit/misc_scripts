import os
import sys
import argparse
import math
from Bio import AlignIO


def parse_and_output_alignment(file_path, start_coord, end_coord, alignment_front_back_buffer):
    # Open the file containing the alignment
    with open(file_path, 'r') as handle:
        # Read the alignment using BioPython's AlignIO
        alignment = AlignIO.read(handle, "stockholm")

        # Extract the reference sequence and remove gaps ('.')
        rf_seq = str(alignment.column_annotations['reference_annotation'])
        rf_seq_0dot = rf_seq.replace('.', '')
        rf_interest_0dot = rf_seq_0dot[start_coord - 1:end_coord]

        # Search for adjusted start coordinate to handle gaps
        rf_seq_len = len(rf_seq)
        start_coord_adj = start_coord
        end_coord_adj = end_coord
        end_coord_adj_found = -1
        while start_coord_adj <= rf_seq_len - (start_coord - end_coord + 1):
            rf_subseg_0dot = rf_seq[start_coord_adj - 1:].replace('.', '')
            start_coord_adj_found = rf_subseg_0dot.find(rf_interest_0dot)
            if start_coord_adj_found != -1 and start_coord_adj_found == 0:
                while end_coord_adj <= rf_seq_len:
                    if rf_seq[start_coord_adj - 1:end_coord_adj].replace('.', '') == rf_interest_0dot:
                        end_coord_adj_found = 1

                        # Extract sequence with buffer
                        seq_spacer = "  ".join([
                            str(alignment[0].seq[start_coord_adj - 1 - alignment_front_back_buffer: start_coord_adj - 1]),
                            str(alignment[0].seq[start_coord_adj - 1:end_coord_adj]),
                            str(alignment[0].seq[end_coord_adj:end_coord_adj + alignment_front_back_buffer])
                        ])

                        # Collect annotations for the aligned sequences
                        column_anno = [seq_spacer]
                        column_anno.append("  ".join([
                            alignment.column_annotations['reference_annotation'][
                            start_coord_adj - 1 - alignment_front_back_buffer: start_coord_adj - 1],
                            alignment.column_annotations['reference_annotation'][start_coord_adj - 1:end_coord_adj],
                            alignment.column_annotations['reference_annotation'][
                            end_coord_adj:end_coord_adj + alignment_front_back_buffer]
                        ]))
                        for key in reversed(list(alignment.column_annotations.keys())):
                            if key != 'reference_annotation':
                                column_anno.append("  ".join([
                                    alignment.column_annotations[key][
                                    start_coord_adj - 1 - alignment_front_back_buffer: start_coord_adj - 1],
                                    alignment.column_annotations[key][start_coord_adj - 1:end_coord_adj],
                                    alignment.column_annotations[key][
                                    end_coord_adj:end_coord_adj + alignment_front_back_buffer]
                                ]))

                        display_len = alignment_front_back_buffer * 2 + (end_coord_adj - start_coord_adj + 1)

                        print_len = 100
                        print_block_num = math.ceil(display_len / print_len)

                        # Print the alignment with appropriate line breaks
                        if print_block_num > 1:
                            for i in range(0, print_block_num):
                                for j in column_anno:
                                    if i == print_block_num - 1:
                                        print(j[100 * i:])
                                    else:
                                        print(j[100 * i:100 * (i + 1)])
                                print()
                        else:
                            for j in column_anno:
                                print(j)
                            print()

                        break
                    else:
                        end_coord_adj += 1
                if end_coord_adj_found != -1:
                    break
            else:
                start_coord_adj = start_coord_adj + start_coord_adj_found
                end_coord_adj = end_coord_adj + start_coord_adj_found


def main(arguments):
    # Define command-line argument parser
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('stk_file', help="Input file", type=str)
    parser.add_argument('coord1', help="reference coordinate 1", type=int)
    parser.add_argument('coord2', help="reference coordinate 2", type=int)

    # Parse command-line arguments
    args = parser.parse_args(arguments)

    # Adjust end coordinate to start coordinate plus 49 if the gap between coordinates is less than 10
    if args.coord2 - args.coord1 < 10:
        args.coord2 = args.coord1 + 49

    alignment_front_back_buffer = 50

    # Call function to parse and output alignment
    parse_and_output_alignment(args.stk_file, args.coord1, args.coord2, alignment_front_back_buffer)


if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))

