import sys
import re

a=sys.argv[1]
b=sys.argv[2]
verbose=""

if a == b: sys.stderr.write( 'INFO:: the infiles are identical will not check if content is identical only validate row and column names.\n' )

while not re.match('^(([Yy]((es)|(ES))?)|([Nn][Oo]?))$',verbose): verbose=raw_input('Show details? (yes/no) ')
if re.match('^[Yy]((es)|(ES))?$',verbose): verbose=True
else: verbose=False

sys.stderr.write( 'Loading input files:\n' )
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

sys.stderr.write( 'Validating row and column names in counts matrix ...\n' )

gene_name_errors = 0
name_in_a_not_in_b = 0
name_in_b_not_in_a = 0
spot_errors = 0
spot_in_a_not_in_b = 0
spot_in_b_not_in_a = 0

if a_row_number != b_row_number:
    if verbose: sys.stderr.write( 'Error:: the number of rows (spots) does not match in the input files.\n' )

_tmp_col_names = dict()
for column_name in a_column_index:
    if column_name not in _tmp_col_names: _tmp_col_names[column_name] = True
    else: sys.stderr.write( 'WARNING:: column (gene) {0} in {1} is present more than once\n'.format(column_name,a) )
    if column_name not in b_column_index:
        if verbose: sys.stderr.write( 'Error:: column (gene) {0} in {1} is not present in {2}\n'.format(column_name,a,b) )
        gene_name_errors += 1
        name_in_a_not_in_b += 1
_tmp_col_names = dict()
for column_name in b_column_index:
    if column_name not in _tmp_col_names: _tmp_col_names[column_name] = True
    else: sys.stderr.write( 'WARNING:: column (gene) {0} in {1} is present more than once\n'.format(column_name,b) )
    if column_name not in a_column_index:
        if verbose: sys.stderr.write( 'Error:: column (gene) {0} in {1} is not present in {2}\n'.format(column_name,b,a) )
        gene_name_errors += 1
        name_in_b_not_in_a += 1


_row_names = dict()
for row_name in a_rows:
    if row_name not in _row_names: _row_names[row_name] = True
    else: sys.stderr.write( 'WARNING:: row (spot) {0} in {1} is present more than once\n'.format(row_name,a) )
    if row_name not in b_rows:
        if verbose: sys.stderr.write( 'Error:: row {0} in {1} is not present in {2}\n'.format(row_name,a,b) )
        spot_errors +=1
        spot_in_a_not_in_b +=1
_row_names = dict()
for row_name in b_rows:
    if row_name not in _row_names: _row_names[row_name] = True
    else: sys.stderr.write( 'WARNING:: row (spot) {0} in {1} is present more than once\n'.format(row_name,b) )
    if row_name not in a_rows:
        if verbose: sys.stderr.write( 'Error:: row {0} in {1} is not present in {2}\n'.format(row_name,b,a) )
        spot_errors +=1
        spot_in_b_not_in_a +=1

resume = "Yes"
if gene_name_errors or spot_errors:
    resume = ""
    sys.stderr.write( 'Error(s) were identified. ' )
    while not re.match('^(([Yy]((es)|(ES))?)|([Nn][Oo]?))$',resume):
        resume=raw_input('continue anyways? (yes/no) ')

if a == b:  sys.exit(0)

if re.match('^[Yy]((es)|(ES))?$',resume):
    sys.stderr.write( 'Searching for value (counts) differences in b table with a as reference ...\n' )
    value_differences = 0
    _tmp_counter = 0
    for row_name, columns in a_rows.iteritems():
    
        sys.stderr.write( str(round(100.0*_tmp_counter/a_row_number,2))+'% Row '+str(row_name)+':\r')
    
        for name, index in a_column_index.iteritems():
    
            a_value = columns[index]
            try:
                b_value = b_rows[row_name][ b_column_index[name] ]
            except KeyError:
                if verbose: sys.stderr.write( '\ncolumn {} notfound in b \n'.format(name) )
                #sys.exit()
    
            if a_value != b_value:
                value_differences += 1
                if verbose: sys.stderr.write( '\nERROR:: the values don\'t match {0} != {1} for {2} at {3} (file {4} vs {5})\n'.format(str( a_value ),str( b_value),name,row_name,a,b) )
                continue
                
            if re.match('^[0-9]+(\.[0-9]+)?$',a_value) or re.match('^[0-9]+(\.[0-9]+)?$',b_value):
                if float(a_value) != float(b_value):
                    value_differences += 1
                    if verbose: sys.stderr.write( '\nERROR:: the values (floats) dont match {0} != {1} for {2} at {3} (file {4} vs {5})\n'.format(str( a_value ),str( b_value),name,row_name,a,b) )
    
        _tmp_counter+=1

print "\n\n:: Summary ::\n"
print "A total of {0} rows (spots) and {1} genes was found in {2}".format(len(a_rows),len(a_column_index),a)
print "A total of {0} rows (spots) and {1} genes was found in {2}".format(len(b_rows),len(b_column_index),b)
print "\n{0} gene name differences and {1} spot coordinate differences was observed between the two input files.".format(gene_name_errors,spot_errors)
print "    {0} gene name(s) and {1} spot coordinate(s) was observed in {2} but not in {3}.".format(name_in_a_not_in_b,spot_in_a_not_in_b,a,b)
print "    {0} gene name(s) and {1} spot coordinate(s) was observed in {2} but not in {3}.".format(name_in_b_not_in_a,spot_in_b_not_in_a,b,a)
if re.match('^[Yy]((es)|(ES))?$',resume):
    print "\n{0} values differed for spots and genes with the same ids in the two counts matrixes.".format(value_differences)

if gene_name_errors == 0 and spot_errors == 0 and value_differences == 0:
    print '\nAll values identical in count tabels.'
else:
    print '\nThe counts matrixes {0} and {1} ARE NOT THE SAME.'.format(a,b)