import sys
import gzip
from estimate_size import estimate_read_count
import time
import itertools
import os
import re
import random
import math

def get_reads(file_handle_1,file_handle_2):
    _buf_1 =[]
    _buf_2 =[]
    for line1,line2 in itertools.izip(file_handle_1,file_handle_2):
        _buf_1.append(line1)
        _buf_2.append(line2)
        if len(_buf_1) == 4:
            yield _buf_1,_buf_2
            _buf_1 =[]
            _buf_2 =[]

def validate(filename, stop_after=10):

    mates_filename = ''
    if   filename.count('_R1_') == 1 and filename.count('_R2_') == 0:  mates_filename = filename.replace('_R1_','_R2_')
    elif filename.count('_R2_') == 1 and filename.count('_R1_') == 0:mates_filename = filename.replace('_R2_','_R1_')
    elif filename.count('_1')==1 and filename.count('_2')==0: mates_filename = filename.replace('_1','_2')
    elif filename.count('_1')==0 and filename.count('_2')==1: mates_filename = filename.replace('_2','_1')
    else:
        mates_filename=raw_input('INPUT:: The input file name is {} the software can\'t automatically detect the file with the mates, please supply the mates file name here: '.format(filename))
    
    assert os.path.exists(mates_filename), 'ERROR:: the file with mates can\'t be found.\n'
    
    print 'INFO:: Will validate as many pairs as possible during {} seconds'.format(stop_after)

    print 'INFO:: Opening files ....'
    if filename.endswith('.gz') and mates_filename.endswith('.gz'):
        open_func = gzip.open
    elif mates_filename.endswith('.'+filename.split('.')[-1]):
        open_func = open
    else:
        print 'ERROR:: mates don\'t hade the same extension.'
        sys.exit(1)
    reads_file = open_func(filename)
    reads_file2 = open_func(mates_filename)
    
    counter_all = 0
    start_time = time.time()
    int_time =0
    checked=0
    print 'INFO:: R1={}'.format(filename)
    print 'INFO:: R2={}'.format(mates_filename)
    print 'INFO:: validating infiles...'
    for number,read_pair in enumerate(get_reads(reads_file,reads_file2)):
        counter_all += 1
        
        if random.randint(0,100) <= 0: # check one read in thousand
            checked+=1
            hed1=read_pair[0][0].rstrip()
            hed2=read_pair[1][0].rstrip()
            seq1=read_pair[0][1].rstrip()
            seq2=read_pair[1][1].rstrip()
            sep1=read_pair[0][2].rstrip()
            sep2=read_pair[1][2].rstrip()
            qul1=read_pair[0][3].rstrip()
            qul2=read_pair[1][3].rstrip()
            assert hed1 == hed2, "\nERROR:: the R1 R2 headers don't match!\n"
            assert re.match("^[AGTCN]+$",seq1), "\nERROR:: non DNA/RNA sequence in R1!\n"
            assert re.match("^[AGTCN]+$",seq2), "\nERROR:: non DNA/RNA sequence in R2!\n"
            assert sep1 == '+', "\nERROR:: R1 qual sequence sepparator is not '+'!\n"
            assert sep2 == '+', "\nERROR:: R1 qual sequence sepparator is not '+'!\n"
            assert re.match("^[\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>\?\@ABCDEFGHIJ]+$",qul1), "\nERROR:: non sanger qual value in R1!\n"
            assert re.match("^[\!\"\#\$\%\&\'\(\)\*\+\,\-\.\/0123456789\:\;\<\=\>\?\@ABCDEFGHIJ]+$",qul2), "\nERROR:: non sanger qual value in R2!\n"
            assert len(seq1) == len(qul1), "\nERROR:: seq and qual length diff in R1!\n"
            assert len(seq2) == len(qul2), "\nERROR:: seq and qual length diff in R2!\n"
        
        #if counter_all%int(round(estimated_read_count/10000.0,0))==0:
        # if counter_all%10000==0:
            # sys.stderr.write(
            #     'INFO:: Progress {}% of input file completed, ETA={}s\r'.format(
            #         'NA',#round(100*float(counter_all)/estimated_read_count,2),
            #         round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)
            #         )
            #     )
        if time.time()-start_time-1 > int_time:
            int_time = int(math.ceil(time.time()-start_time))-1
            sys.stderr.write(
                'INFO:: Seconds since start {}, {} read pairs processed, {} random read pairs validated\r'.format(
                    int_time,
                    counter_all,
                    checked
                    )
                )

        if time.time()-start_time >= stop_after: break
    print ''
    
    print 'INFO:: Parsed a total of {} reads from infile.'.format(counter_all)
    
    print 'INFO:: closing files...'
    reads_file.close()
    reads_file2.close()
    
    print 'INFO:: completed.'

if __name__ == '__main__':
    validate(sys.argv[1])