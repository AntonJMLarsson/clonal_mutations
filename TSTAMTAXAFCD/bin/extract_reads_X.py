import argparse
import pysam

def main():
    parser = argparse.ArgumentParser(description='Subset bam file by contig and only include exon-mapping reads')
    parser.add_argument('-i','--input', type=str, help='Input .bam file')
    parser.add_argument('-o','--output', type=str, help='Output .bam file')
    parser.add_argument('-s','--suffix', type=str, default='', help='Modify cell barcode with suffix')
    parser.add_argument('--contig', default='X', metavar='contig', type=str, help='Restrict to contig')
    parser.add_argument('--cells', default=None, type=str, help='Cell barcodes')

    args = parser.parse_args()
    bam_in = pysam.AlignmentFile(args.input, 'rb')
    bam_out = pysam.AlignmentFile(args.output, 'wb', template=bam_in)

    cell_set = set([line.rstrip() for line in open(args.cells)])

    for read in bam_in.fetch(args.contig):
        if (read.has_tag('GE') or read.has_tag('GI')) and read.get_tag('BC') in cell_set:
            if args.suffix != '':
                cell = read.get_tag('BC')
                cell = '{}_{}'.format(cell, args.suffix)
                read.set_tag('BC', cell)
                bam_out.write(read)
            else:
                bam_out.write(read)
    bam_in.close()
    bam_out.close()

if __name__ == '__main__':
    main()
