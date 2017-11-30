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
    # Reduce total read count by 10% as the estimate function usually overestimates this means we might alway skip the last part of the file
    # estimated_read_count = int(round(estimated_read_count*0.9,0))
    
    print 'selecting reads...'
    selected = {}
    
    mode='select'
    if subset_size >= estimated_read_count:
        print 'ERROR:: subset is larger than estimade file size exiting'
        sys.exit()
    # Mode exclude is not fully implemented
    # elif subset_size >= estimated_read_count*0.9:
    #     mode='exclude'
    #     print 'WARNING:: Subset is greater than 90% of input reads forcing mode=exclude size of the subset written can\'t be guaranteed'
    #     subset_size = estimated_read_count-subset_size

    rand_min = 0
    rand_max = estimated_read_count
    start_time = time.time()
    while len(selected) < subset_size:
        # if selected and len(selected)%int(round(subset_size/10.0,0))==0:
        #     sorted_keys = sorted(selected.keys())
        #     print ''
        #     if sorted_keys[0] == rand_min:
        #         i = 0
        #         while sorted_keys[i]+1 == sorted_keys[i+1]: i+=1
        #         print 'pos',i,'val',sorted_keys[i],'next',sorted_keys[i]+1, '!=', sorted_keys[i+1]
        #         rand_min = sorted_keys[i]
        #     if sorted_keys[-1] == rand_max:
        #         i=-1
        #         while sorted_keys[i]-1 == sorted_keys[i-1]: i-=1
        #         print 'pos',i,'val',sorted_keys[i],'next',sorted_keys[i]-1, '!=', sorted_keys[i-1]
        #         rand_max=sorted_keys[i]
        #     print 'getting random integers between {} and {}'.format(rand_min,rand_max)
        key = random.randint(rand_min,rand_max)
        selected[key] = False
        if len(selected)%int(round(subset_size/1000,0))==0 or len(selected)==1 or time.time()-last_time>1:
            last_time = time.time()
            sys.stderr.write(
                'Progress selected {}% of subset size, ETA={}s\r'.format(
                            round(100*float(len(selected))/(subset_size),2),
                            round((time.time()-start_time)*subset_size/float(len(selected))-(time.time()-start_time), 2)
                        )
                    )
    print ''
    
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
    print 'copying selected reads from infile to outfile...'
    for number,read in enumerate(get_reads(reads_file)):
        counter_all += 1
        if counter_all%int(round(estimated_read_count/1000,0))==0: sys.stderr.write('Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(round(100*float(counter_all)/estimated_read_count,1), round(100*float(counter_mode)/subset_size,1), round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)))
        if mode=='select' and number in selected:
            outfile.write( ''.join(read) )
            selected[number] = True
            counter_mode+=1
            counter_written+=1
            if counter_mode%int(round(subset_size/1000,0))==0: sys.stderr.write('Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(round(100*float(counter_all)/estimated_read_count,1), round(100*float(counter_mode)/subset_size,1), round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)))
        elif mode=='select' and number not in selected:
            if random.randint(0,10)==1: extra_reads.append(read)
        elif mode=='exclude' and number in selected:
            selected[number] = True
            counter_mode+=1
            if counter_mode%int(round(subset_size/1000,0))==0: sys.stderr.write('Progress {}% of input file and {}% of subset completed, ETA={}s\r'.format(round(100*float(counter_all)/estimated_read_count,1), round(100*float(counter_mode)/subset_size,1), round((time.time()-start_time)*estimated_read_count/float(counter_all)-(time.time()-start_time),2)))
        elif mode=='exclude' and number not in selected:
            outfile.write( ''.join(read) )
            counter_written+=1
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
    print sum([1 for i in selected.values() if i]),'of the',len(selected),'selection/exlution-dict reads were processed from the outfile.'
    
    print 'closing files...'
    reads_file.close()
    outfile.close()
    
    print 'completed'

if __name__ == '__main__':
    random_subset(sys.argv[1])