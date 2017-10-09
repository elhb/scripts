import sys
import re

a=sys.argv[1]
b=sys.argv[2]

sys.stderr.write( 'loading a ...\n' )
a_column_index = dict()
a_rows = dict()
with open(a) as input_file:
    a_row_number = 0
    for row in input_file:
        columns=row.rstrip().split('\t')
        if a_row_number == 0:
            for index in range(len(columns)):
                a_column_index[columns[index]] = index
        else:
            a_rows[columns[0]] = columns
        a_row_number += 1

sys.stderr.write( 'loading b ...\n' )
b_column_index = dict()
b_rows = dict()
with open(b) as input_file:
    b_row_number = 0
    for row in input_file:
        columns=row.rstrip().split('\t')
        if b_row_number == 0:
            for index in range(len(columns)):
                b_column_index[columns[index]] = index
        else:
            b_rows[columns[0]] = columns
        b_row_number += 1

sys.stderr.write( 'Validating row and column names\n' )
assert a_row_number == b_row_number, 'Error:: the number of rows (spots) does not match'
for row_name in a_rows: assert row_name in b_rows, 'Error:: row {0} in {1} is not present in {2}'.format(row_name,a,b)
for row_name in b_rows: assert row_name in a_rows, 'Error:: row {0} in {1} is not present in {2}'.format(row_name,b,a)
for column_name in a_column_index: assert column_name in b_column_index,  'Error:: column (gene) {0} in {1} is not present in {2}'.format(column_name,a,b)
for column_name in b_column_index: assert column_name in a_column_index,  'Error:: column (gene) {0} in {1} is not present in {2}'.format(column_name,b,a)

sys.stderr.write( 'Searching for differences in b table with a as reference ...\n' )

_tmp_counter = 0
for row_name, columns in a_rows.iteritems():
    #sys.stderr.write( str(round(100.0*_tmp_counter/a_row_number,2))+'% Row '+str(row_name)+':\n')
    sys.stderr.write( str(round(100.0*_tmp_counter/a_row_number,2))+'% Row '+str(row_name)+':\r')
    for name, index in a_column_index.iteritems():
        #sys.stderr.write( ' column '+name)
        a_value = columns[index]
        #sys.stderr.write( str(a_value) )
        try:
            b_value = b_rows[row_name][ b_column_index[name] ]
        except KeyError:
            sys.stderr.write( '\ncolumn {} notfound in b \n'.format(name) )
            sys.exit()
        #sys.stderr.write( str(b_value) )
        try:
            assert a_value == b_value
            if re.match('^[0-9]+(\.[0-9]+)?$',a_value) or re.match('^[0-9]+(\.[0-9]+)?$',b_value):
                assert float(a_value) == float(b_value)
        except:
            sys.stderr.write( '\nERROR:: {0} != {1} for {2} at {3} (file {4} vs {5})\n'.format(str( a_value ),str( b_value),name,row_name,a,b) )
            sys.exit()
        #sys.stderr.write( ' OK                                              \r')
    _tmp_counter+=1


print 'All values identical in count tabels.'