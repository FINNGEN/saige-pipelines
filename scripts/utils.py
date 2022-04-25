import os,subprocess,sys,csv,gzip,logging,mmap,re
from tempfile import NamedTemporaryFile
from functools import partial

def return_columns(l,columns):
    '''
    Returns all columns, or rather the elements, provided the columns
    '''
    if columns == 'all':
        return l
    elif type(columns) == int:
        return l[columns]
    elif type(columns) == list:
        return list(map(l.__getitem__,columns))

    
def pretty_print(string,l = 30):
    l = l-int(len(string)/2)
    print('-'*l + '> ' + string + ' <' + '-'*l)

def tmp_bash(cmd,check = False):
    

    scriptFile = NamedTemporaryFile(delete=True)
    with open(scriptFile.name, 'w') as f:
        f.write("#!/bin/bash\n")
        f.write(cmd + "\n")

    os.chmod(scriptFile.name,0o777)
    scriptFile.file.close()

    if check:
        subprocess.check_call(scriptFile.name)
    else:
        subprocess.call(scriptFile.name,stderr = subprocess.DEVNULL)
        
def make_sure_path_exists(path):
    import errno
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise
        
def file_exists(fname):
    '''
    Function to pass to type in argparse
    '''
    if os.path.isfile(fname):
        return str(fname)
    else:
        print(fname + ' does not exist')
        sys.exit(1)


def basic_iterator(f,separator = None,skiprows = 0,count = False,columns = 'all'):
    '''
    Function that iterates through a file and returns each line as a list with separator being used to split.
    '''

    # converts column names to indexes
    if type(columns) == list:
        if all(isinstance(n, str) for n in columns):
            header = return_header(f)
            columns = [header.index(elem) for elem in columns]
   
    open_func = return_open_func(f)
    if not separator:separator = identify_separator(f)
    i = open_func(f)
    for x in range(skiprows):next(i)

    if count is False:
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)          
            yield line
    else:
        row = 0
        for line in i:
            line =line.strip().split(separator)
            line = return_columns(line,columns)
            row += 1   
            yield row,line


def return_open_func(f):
    '''
    Detects file extension and return proper open_func
    '''
   
    file_path,file_root,file_extension = get_path_info(f)

    if 'bgz' in file_extension:
        #print('gzip.open with rb mode')
        open_func = partial(gzip.open, mode = 'rb')
    
    elif 'gz' in file_extension:
        #print('gzip.open with rt mode')
        open_func = partial(gzip.open, mode = 'rt')

    else:
        #print('regular open')
        open_func = open      
    return open_func



def identify_separator(f):
    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    sniffer = csv.Sniffer()
    dialect = sniffer.sniff(header)
    return dialect.delimiter


def get_path_info(path):
    file_path = os.path.dirname(path)
    basename = os.path.basename(path)
    file_root, file_extension = os.path.splitext(basename)
    return file_path,file_root,file_extension


def return_header(f):

    open_func = return_open_func(f)
    with open_func(f) as i:header = i.readline().strip()
    delimiter = identify_separator(f)
    header = header.split(delimiter)
    return header

log_levels = {
    'critical': logging.CRITICAL,
    'error': logging.ERROR,
    'warn': logging.WARNING,
    'warning': logging.WARNING,
    'info': logging.INFO,
    'debug': logging.DEBUG
}

def extract_int_from_string(s):
    return re.search(r'\d+', s).group()


def progressBar(value, endvalue, bar_length=20):
    '''
    Writes progress bar, given value (eg.current row) and endvalue(eg. total number of rows)
    '''

    percent = float(value) / endvalue
    arrow = '-' * int(round(percent * bar_length)-1) + '>'
    spaces = ' ' * (bar_length - len(arrow))

    sys.stdout.write("\rPercent: [{0}] {1}%".format(arrow + spaces, int(round(percent * 100))))
    sys.stdout.flush()
    
def check_region(locus,region):
    if ':' in locus and  '-' in locus:
        locus,region = region,locus
    return locus,f" --range {region} "
