#!/usr/bin/env python3
import pandas as pd
import argparse, urllib




parser = argparse.ArgumentParser(description="""

    Input
    - List of Fbn dmel gene names to query FlyAtlas2 with.

    Returns
    - table with genes in rows and FPKMs for all tissues, sexes, and stages in columns.
    - error file with names of genes that were not added to above table
    - raw atlas tables text file - enabled URL-free searches in the future as described below.

    
    - FUTURE / TODO:
        - optional return:
        - returns gzipped tarball with all individual gene tables (of given gene names) that can be used as input for future queries rather accessing the URL.
        - Or as single gzipped text file where entries are separated w/ recognizable lines.
        - in both cases, for querying, one would need to collect all gene names first,
            - then do a single pass over the tarball or text file.
    """, formatter_class= argparse.RawTextHelpFormatter)



parser.add_argument('-f','--genefile',
                   type=str, default=False,
                   help='''Path to input file of Fbn Dmel gene names - 1 per line, single column. You can use this alone or with -c/--genes.''')
parser.add_argument('-c','--genes',
                   type=str, default=False,
                   help='''Give comma=sep list of Fbn gene names at commandline. You can use this alone or with -f/--genefile.''')
parser.add_argument('-o','--outprefix',
                   type=str, default='myFlyAtlas2Queries',
                   help='''Prefix for output files.
Default = myFlyAtlas2Queries.
Careful : Using default for multiple applications of this script in the same directory will result in over-writing each previous output. ''')


args = parser.parse_args()



## FUNCTIONS and VARS
def mackerel(g):
    # mackerel just cute way of saying "make url"
    return 'http://flyatlas.gla.ac.uk/FA2Direct/index.html?fbgn=' + g + '&tableOut=gene'
    
head = ['ID', 'anno','sym','name']
tissues = ['Whole body', 'Head', 'Eye', 'Brain / CNS', 'Thoracicoabdominal ganglion', 'Crop',
          'Midgut', 'Hindgut', 'Malpighian Tubules', 'Fat body', 'Salivary gland',
          'Heart', 'Trachea', 'Ovary', 'Virgin Spermatheca', 'Mated Spermatheca', 'Testis',
          'Accessory glands', 'Carcass', 'Rectal pad']
fpkm_idx = {'male':1, 'female':4, 'larval':9}

def initialize():
    l = []
    l += head
    for t in tissues:
        for stg in fpkm_idx.keys():
            key = stg + '_' + '_'.join([e for e in t.split() if e not in ['/']])
            l.append(key)
    return pd.DataFrame(data=None, columns=l)

def wrangle(g):
    url = mackerel(g)
    con = urllib.request.urlopen(url)
    a = con.read()
    stringtable = a.decode('utf-8')
    data = [e.split('\t') for e in stringtable.strip().split('\n')]
    return data, stringtable

def get_body_idx(data):
    for i in range(len(data)):
        if data[i][0] == 'Whole body':
            return i
    print("Expectation Error: Whole body not found in table....")
    return

        
def process(data):
    bi = get_body_idx(data)
    header = data[:bi]
    body = {e[0]:e for e in data[bi:]}
    try:
        assert set(tissues) == set(body.keys())
    except AssertionError:
        print(tissues)
        print (body.keys())
    d = {}
    d['ID'] = header[0][-1]
    d['anno'] = header[1][-1]
    d['sym'] = header[2][-1]
    d['name'] = header[3][-1]
    for t in tissues:
        for stg in fpkm_idx.keys():
            key = stg + '_' + '_'.join([e for e in t.split() if e not in ['/']])
            d[key] = body[t][fpkm_idx[stg]]
    return d

def get_gene_list_f(fh):
    '''fh goes to file of flybase gene names (1 per line) - assumes single column'''
    return [e.strip() for e in open(fh).readlines()]

def get_gene_list_c(s):
    '''fh goes to file of flybase gene names (1 per line) - assumes single column'''
    return [e.strip() for e in s.strip().split(',')]


def run_pipeline(args):
    '''fh goes to file of flybase gene names (1 per line) - assumes first column'''
    # Get all gene names provided; make sure they're unique
    genes = []
    if args.genefile:
        genes += get_gene_list_f(args.genefile)
    if args.genes:
        genes += get_gene_list_c(args.genes)
    genes = list(set(genes))

    # Init dataframe
    df = initialize()
    singleText = ''
    count = -1
    errors = []
    for g in genes:
        count += 1
        if count%100 == 0:
            print("Current gene number is:", count)
            print("Current gene is:", g)
        try:
            data, stringtable = wrangle(g)
            singleText += 'NewEntry\n' + stringtable + '\n'*2
            d = process(data)
            df = df.append(d, ignore_index=True)
        except:
            print("Error on:", g)
            errors.append(g)
        
    return df, errors, singleText


def run(args):
    df, errors, singleText = run_pipeline(args)

    # write
    # Dataframe
    df.to_csv(args.outprefix + '.txt', sep="\t", index=False)
    # Errors
    with open(args.outprefix + '-errors.txt','w') as f:
        for e in errors:
            f.write(e+'\n')
    # Single Text File
    with open(args.outprefix + '-rawAtlasTables.txt','w') as g:
        g.write(singleText)





############ EXECUTE

run(args)

