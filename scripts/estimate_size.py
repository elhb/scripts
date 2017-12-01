import sys
import gzip
import os
import time
import random

def estimate_read_count(filename, gzipcompressionlevel=1):
    """
    Function used to get a rough estimate o readcount in fastq files.
    It reads a subset of 100 or 200 MB (depending on compression)
    and counts the number of lines in the subset. The line count, subset size
    and the size of full file is then used to calulate an estimated read count.
    This function often seems to overestimate the read count by up to 10%.
    
    Requiered Input Parameters:
       filename :: name of the input file
    Optional:
       gzipcompressionlevel :: Gzip compression level to use when comparing
                               file sizes. Default is 1 as this seems to be
                               standard for files straigh of the MiniSeq.
    Returns:
       The estimated total readcount as an integer
    """

    if filename.endswith('.gz'):
        open_func = gzip.open
        subset_size = 100*1024**2
    else:
        open_func = open
        subset_size = 200*1024**2
    reads_file = open_func(filename)
    
    data = []
    counter=0
    data=reads_file.read(subset_size)
    
    counter = data.count('\n')
    
    tmp_file_name='DELETEME_IM_A_TEMPORARY_FILE.'+str(''.join([str(random.randint(0,9)) for i in xrange(20)]))
    if filename.endswith('.gz'):
        tmp_file = gzip.open(tmp_file_name,'w', compresslevel=gzipcompressionlevel)
    else:
        tmp_file = open(tmp_file_name,'w')
    
    tmp_file.write( data )
    tmp_file.close()
    fraction=float(os.path.getsize(filename))/os.path.getsize(tmp_file_name)
    os.remove(tmp_file_name)
    
    estimate = int(round(float(counter)*fraction/4.0,0))
    if __name__ == '__main__':
        print 'Total read count estimated to ~',int(round(float(counter)*fraction/4.0,0)),'reads.'
        print 'The estimate was based on a subset of',counter,'lines.'
    
    return estimate

if __name__ == '__main__':
    estimate_read_count(sys.argv[1])