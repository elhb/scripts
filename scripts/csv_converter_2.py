import sys
import csv

def converter(csv_filename):
    """
    Small converter script to be used to convert the output of dstat
    to a format that is easier for the eye and excel or google sheets
    Usage:
        python csv_converter_2.py <infile.csv>
    """
    
    infile = csv.reader(open(csv_filename))

    # Example of infeile format:
    # "Dstat 0.7.3 CSV output"
    # "Author:","Dag Wieers <dag@wieers.com>",,,,"URL:","http://dag.wieers.com/home-made/dstat/"
    # "Host:","stlinuxsrv01",,,,"User:","eb"
    # "Cmdline:","dstat -Tmcd --top-cpu --noheader --output results_171011_12_36/MOB_REPLICATE_5_pipe_filerInput/sys_stat.MOB_REPLICATE_5_pipe_filerInput.csv",,,,"Date:","11 Oct 20
    # "epoch","memory usage",,,,"total cpu usage",,,,,"dsk/total",,"most expensive"
    # "epoch","used","free","buff","cach","usr","sys","idl","wai","stl","read","writ","cpu process"

    # trash the first 5 lines
    for i in range(5): infile.next()
    
    #get the header
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
        try: process = line[header_dict['cpu process']].split(' / ')[0]
        except:
            sys.stderr.write(
                "WARNING:: could not identify the cpu process for line {}\n".format(",".join(line))
                )
            process= 'OTHER'
        data.append(
                (
                    float(line[header_dict['epoch']]),
                    round(float(line[header_dict['used']])/1024**2,2),
                    round(float(line[header_dict['usr']]),2),
                    float(line[header_dict['read']]),
                    float(line[header_dict['writ']]),
                    110 if process == 'STAR' else 200,
                    110 if process == 'taggd_demultipl' else 200,
                    110 if process == 'st_pipeline_run' else 200
                )
            )
        line_counter+=1
    
    for line in data:
        _1,_2,_3,_4,_5,_6,_7,_8 = line
        _1 = str(_1-min([_tmp[0] for _tmp in data])).replace('.',',')
        _2 = str(_2).replace('.',',')
        _3 = str(_3).replace('.',',')
        _4 = str(round(_4/max([_tmp[3] for _tmp in data])*100,2)).replace('.',',')
        _5 = str(round(_5/max([_tmp[4] for _tmp in data])*100,2)).replace('.',',')
        _6 = str(_6).replace('.',',')
        _7 = str(_7).replace('.',',')
        _8 = str(_8).replace('.',',')
        print "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(_1,_2,_3,_4,_5,_6,_7,_8)

if __name__ == '__main__':
    converter(sys.argv[1])