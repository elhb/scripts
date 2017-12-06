import sys
import random
import gzip
from estimate_size import estimate_read_count
import time
import itertools
import os
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

def random_subset(filename,subset_size=int(round(3e5,0))):

    mates_filename = ''
    if   filename.count('_R1_') == 1 and filename.count('_R2_') == 0:  mates_filename = filename.replace('_R1_','_R2_')
    elif filename.count('_R2_') == 1 and filename.count('_R1_') == 0:mates_filename = filename.replace('_R2_','_R1_')
    elif filename.count('_1')==1 and filename.count('_2')==0: mates_filename = filename.replace('_1','_2')
    elif filename.count('_1')==0 and filename.count('_2')==1: mates_filename = filename.replace('_2','_1')
    else:
        mates_filename=raw_input('INPUT:: The input file name is {} the software can\'t automatically detect the file with the mates, please supply the mates file name here: '.format(filename))
    
    assert os.path.exists(mates_filename), 'ERROR:: the file with mates can\'t be found.\n'
    
    original_subset_size = subset_size
    outfilename = filename+'.random_subset.{}_reads.fq'.format(original_subset_size)+{True:'.gz',False:''}[filename.endswith('.gz')]
    mates_outfilename = mates_filename+'.random_subset.{}_reads.fq'.format(original_subset_size)+{True:'.gz',False:''}[mates_filename.endswith('.gz')]
    
    print 'INFO:: Estimating read count...'
    if os.path.getsize(filename) > os.path.getsize(mates_filename): estimated_read_count = estimate_read_count(filename)
    else: estimated_read_count = estimate_read_count(mates_filename)
    print 'INFO:: Read count estimated to {}.'.format(estimated_read_count)
    
    mode='select'
    if subset_size >= estimated_read_count:
        print 'ERROR:: subset is larger than estimade file size exiting'
        sys.exit()

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
    print 'INFO:: output file name', outfilename
    print 'INFO:: mates output file name', mates_outfilename
    if filename.endswith('.gz'):
        outfile = gzip.open(outfilename,'w')
        outfile2 = gzip.open(mates_outfilename,'w')
    else:
        outfile = open(outfilename,'w')
        outfile2 = open(mates_outfilename,'w')
    
    counter_all = 0
    counter_mode = 0
    counter_written = 0

    start_time = time.time()
    extra_reads = []
    probablility_base=estimated_read_count
    print 'INFO:: copying selected reads from infile to outfile...'
    for number,read_pair in enumerate(get_reads(reads_file,reads_file2)):
        counter_all += 1
        if counter_all%int(round(estimated_read_count/10000.0,0))==0:
            sys.stderr.write(
                'INFO:: Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(
                    round(100*float(counter_all)/estimated_read_count,2),
                    round(100*float(counter_mode)/subset_size,2),
                    round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)
                    )
                )
        if random.randint(0,int(probablility_base)) <= probablility_base*float(subset_size)*1.05/estimated_read_count and counter_written!=subset_size:
            outfile.write( ''.join(read_pair[0]) )
            outfile2.write( ''.join(read_pair[1]) )
            counter_mode+=1
            counter_written+=1
            if counter_mode%int(round(subset_size/10000.0,0))==0:
                sys.stderr.write(
                    'INFO:: Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(
                        round(100*float(counter_all)/estimated_read_count,2),
                        round(100*float(counter_mode)/subset_size,2),
                        round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)
                        )
                    )
        else:
            # keep a buffer approx 10% the size of the subset to have just in case after all has been written
            if random.randint(0,int(probablility_base)) < 0.1*probablility_base*float(subset_size)/estimated_read_count:
                extra_reads.append(read_pair)
    print ''
    
    print 'INFO:: Parsed a total of {} reads from infile.'.format(counter_all)
    print 'INFO:: Estimated size was {}'.format(estimated_read_count)
    print 'INFO:: Estimation was {}% of actual size.'.format(round(100*float(estimated_read_count)/counter_all,2))
    print 'INFO:: {} of the {} reads were written to the outfile ({}%).'.format(counter_written, subset_size, round(100*float(counter_written)/subset_size,2))
    if counter_written < original_subset_size:
        print 'WARNING:: to few reads were written will use a random selection of already buffered extra reads to fill the subset ...'
        extra_to_write=original_subset_size-counter_written
        if len(extra_reads) < extra_to_write:
            print 'ERROR:: wont be able to fill subset to few reads were saved in buffer (there are {} we need {}).'.fomat(len(extra_reads), extra_to_write)
        else:
            start_time = time.time()
            random.shuffle(extra_reads)
            extra_written = 0
            for read_pair in extra_reads:
                outfile.write( ''.join(read_pair[0]) )
                outfile2.write(''.join(read_pair[1]) )
                counter_written+=1
                extra_written+=1
                if extra_written%int(math.ceil(extra_to_write/10000.0))==0:
                    sys.stderr.write(
                        'INFO:: Progress {}% of extra read pair writes, ETA={}s\r'.format(
                            round(100*float(extra_written)/extra_to_write,2),
                            round((time.time()-start_time)*extra_to_write/float(extra_written)-(time.time()-start_time),2)
                            )
                        )
                if extra_written == extra_to_write: break
        print 'INFO:: {} of the {} reads were written to the outfile.'.format(counter_written,original_subset_size)
    
    print 'INFO:: closing files...'
    reads_file.close()
    outfile.close()
    reads_file2.close()
    outfile2.close()
    
    print 'INFO:: completed.'

if __name__ == '__main__':
    random_subset(sys.argv[1])