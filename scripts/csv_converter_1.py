"""
Small converter script to be used to convert the output of glances
to a format that is easier for the eye and excel or google sheets
Usage:
    python csv_converter_1.py <infile.csv>
"""

import sys
import csv

infile = csv.reader(open(sys.argv[1]))
header_list = infile.next()
header_dict = { header_list[i]:i for i in range(len(header_list)) }

print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(
    'TIME',
    'Memory (MB)',
    'CPU (%)',
    'Read (%)',
    'Write (%)',
    'STAR',
    'TAGGD',
    'st_pipeline',
)
line_counter = 0
data = []
for line in infile:
    try:read = float(line[header_dict['diskio_disk0_read_bytes' ]])
    except ValueError: read = 0
    write= float(line[header_dict['diskio_disk0_write_bytes']])
    data.append((
        line_counter,
        round(float(line[header_dict['mem_used']])/1024**2,2),
        round(sum([ float(line[header_dict[column_name]]) for column_name in ['percpu_0_total','percpu_1_total','percpu_2_total','percpu_3_total'] ])/4.0,2),
        read,
        write,
        200,
        200,
        200))
    line_counter+=1

for line in data:
    _1,_2,_3,_4,_5,_6,_7,_8 = line
    _1 = str(_1).replace('.',',')
    _2 = str(_2).replace('.',',')
    _3 = str(_3).replace('.',',')
    _4 = str(round(_4/max([_tmp[3] for _tmp in data])*100,2)).replace('.',',')
    _5 = str(round(_5/max([_tmp[4] for _tmp in data])*100,2)).replace('.',',')
    _6 = str(_6).replace('.',',')
    _7 = str(_7).replace('.',',')
    _8 = str(_8).replace('.',',')
    print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(_1,_2,_3,_4,_5,_6,_7,_8)
