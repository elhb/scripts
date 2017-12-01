import sys
import random
import gzip
from estimate_size import estimate_read_count
import time

def get_reads(file_handle):
    _buf =[]
    for line in file_handle:
        _buf.append(line)
        if len(_buf) == 4:
            yield _buf
            _buf = []

def random_subset(filename,subset_size=int(round(3e6,0))):

    original_subset_size = subset_size
    outfilename = filename+'.random_subset.{}_reads.fq'.format(original_subset_size)+{True:'.gz',False:''}[filename.endswith('.gz')]
    
    print 'estimating read count...'
    estimated_read_count = estimate_read_count(filename)
    print 'read count estimated to',estimated_read_count
    
    mode='select'
    if subset_size >= estimated_read_count:
        print 'ERROR:: subset is larger than estimade file size exiting'
        sys.exit()

    print 'opening files ....'
    if filename.endswith('.gz'):
        open_func = gzip.open
    else:
        open_func = open
    reads_file = open_func(filename)
    print 'output file name', outfilename
    if filename.endswith('.gz'):
        outfile = gzip.open(outfilename,'w')
    else:
        outfile = open(outfilename,'w')
    
    counter_all = 0
    counter_mode = 0
    counter_written = 0

    start_time = time.time()
    extra_reads = []
    probablility_base=1e9
    print 'copying selected reads from infile to outfile...'
    for number,read in enumerate(get_reads(reads_file)):
        counter_all += 1
        if counter_all%int(round(estimated_read_count/10000,0))==0:
            sys.stderr.write(
                'Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(
                    round(100*float(counter_all)/estimated_read_count,2),
                    round(100*float(counter_mode)/subset_size,2),
                    round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)
                    )
                )
        if random.randint(0,int(probablility_base)) < probablility_base*float(subset_size)/estimated_read_count and counter_written!=subset_size:
            outfile.write( ''.join(read) )
            counter_mode+=1
            counter_written+=1
            if counter_mode%int(round(subset_size/10000,0))==0:
                sys.stderr.write(
                    'Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(
                        round(100*float(counter_all)/estimated_read_count,2),
                        round(100*float(counter_mode)/subset_size,2),
                        round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)
                        )
                    )
        else:
            if random.randint(0,10)==1: extra_reads.append(read)
    print ''
    
    print 'Parsed a total of',counter_all,'reads from infile.'
    print 'Estimated size was',estimated_read_count
    print 'Estimation was',round(100*float(estimated_read_count)/counter_all,2),'% of actual size.'
    print counter_written,'of the',original_subset_size,' reads were written to the outfile.'
    if counter_written < original_subset_size:
        print 'to few reads were written will use some random extra reads to fill the subset ...'
        extra_to_write=original_subset_size-counter_written
        if len(extra_reads) < extra_to_write: print 'ERROR:: wont be able to fill subset to few reads were saved in buffer (there are',len(extra_reads),'we need',extra_to_write,').'
        else:
            _dict = {i:read for i,read in enumerate(extra_reads)}
            while _dict.values().count(False) != extra_to_write:
                random_read_number = random.choice([i for i,read in _dict.iteritems() if read])
                outfile.write( ''.join(_dict[random_read_number]) )
                _dict[random_read_number] = False
                counter_written+=1
        print counter_written,'of the',original_subset_size,' reads were written to the outfile.'
    
    print 'closing files...'
    reads_file.close()
    outfile.close()
    
    print 'completed'

if __name__ == '__main__':
    random_subset(sys.argv[1])