#!/usr/bin/python
import argparse
import sys

def getIntervals( data ):
    data=data.strip()
    linetab = data.split("\t")
    number_of_columns = len(linetab)
    try:
        chr1, pos1, strand1, chr2, pos2, strand2 = linetab[1:7]
    except ValueError:
        print "Error - unexpected format\n", data.strip()
        
    return {'chr1':chr1, 'pos1':pos1, 'strand1':strand1, 'chr2':chr2, 'pos2':pos2, 'strand2':strand2}


def isOrdered( chr1, pos1, chr2, pos2, chromosomes_list):
    if chromosomes_list.index(chr1) < chromosomes_list.index(chr2):
        return True
    elif chr1 == chr2 and int(pos1) < int(pos2):
        return True
    else:
        return False

def isEqual( chr1, pos1, chr2, pos2):
    if chr1==chr2 and pos1==pos2:
        return True
    else:
        return False

def isDuplicated (data1, data2):
    if data1 == None or data2 == None:
        return False

    if isEqual(data1['chr1'], data1['pos1'], data2['chr1'], data2['pos1']) and isEqual(data1['chr2'], data1['pos2'], data2['chr2'], data2['pos2']):
        return True
    else:
        return False

def  compareDataList( data_list, chromosomes_list ):
    
    #print "--------compare------------"
    #print data_list
    #print "--------------------------"

    min_index=None

    ## Check list of available buffers
    allindex=[]
    for i in range(0,len(data_list)):
        if data_list[i] != None:
            allindex.append(i)
        
    if len(allindex)>0:
        ## init with first data
        min_index=allindex[0]
        min_data=data_list[allindex[0]]
        ## compare with other data
        for i in allindex[1:len(allindex)]:
            cur_data=data_list[i]
            ## Compare first mate
            if not isOrdered(min_data['chr1'], min_data['pos1'], cur_data['chr1'], cur_data['pos1'], chromosomes_list):
                min_data=cur_data
                min_index=i
                ## Compare second mate
                if isEqual(min_data['chr1'], min_data['pos1'], cur_data['chr1'], cur_data['pos1']) and not isOrdered(min_data['chr2'], min_data['pos2'], cur_data['chr2'], cur_data['pos2'], chromosomes_list):
                    min_data=cur_data
                    min_index=i
    return min_index


################################################################################
##
##  __MAIN__
##
################################################################################

## Get Args
parser = argparse.ArgumentParser(description='Merge valid pairs Hi-C product from multiple samples.')
parser.add_argument('--rmdup', dest='rmdup', action='store_const', const=True, default=False, help='Remove duplicates')
parser.add_argument('filelist', metavar='validPairs', nargs='+', help='List of Hi-C valid pairs files. These files are expected to be sorted by chromosome and position')

chromosomes_list=['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chrX', 'chrY']
args = parser.parse_args()
rmdup = args.rmdup
nbfiles = len(args.filelist)

## Open handler on valid pairs files
files_handlers = []
line_list = []
data_list = []
count_dup = 0

for cur_file in args.filelist:
    cur_pos=len(files_handlers)
    ##print cur_pos
    cur_handler=open(cur_file, 'r')
    files_handlers.append(cur_handler)
    
    ## Init data list with first line
    firstline = cur_handler.readline().strip()
    lastposition = cur_handler.tell()
    secondline = cur_handler.readline().strip()

    ## Remove duplicates from the same file
    while (isDuplicated(getIntervals(firstline),getIntervals(secondline))):
        lastposition = cur_handler.tell()
        secondline=cur_handler.readline().strip()

    cur_handler.seek(lastposition)
    line_list.append(firstline)
    data_list.append(getIntervals(line_list[cur_pos]))


index=compareDataList(data_list, chromosomes_list)
prec_data=None
while (index != None):
    ## Get sorted line
    ## Print output and remove duplicates

    if not rmdup:
        print line_list[index]
        prec_data=data_list[index]
    else:
        if not isDuplicated(prec_data, data_list[index]):
            print line_list[index]
            prec_data=data_list[index]
        else:
            count_dup+=1

    ## Read a new line for this handler
    newline=files_handlers[index].readline().strip()
    ## Remove duplicates from the same file
    while (newline == line_list[index]):
        newline=files_handlers[index].readline().strip()
    line_list[index]=newline

    if line_list[index]:
        data_list[index]=getIntervals(line_list[index])
    else:
        files_handlers[index].close()
        data_list[index]=None
        line_list[index]=None
            
    index=compareDataList(data_list, chromosomes_list)

print >> sys.stderr, "rm_duplicates\t", count_dup
