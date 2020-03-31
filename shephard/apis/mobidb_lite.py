
## Mobidb-lite input specification
## 
## each line
##
## uniprotID | start | end | sequence
##



def read_in_mobidblite_file(filename):

    with open(filename,'r') as fh:
        content = fh.readlines()

    for line in content:
        if len(line) > 0:
            sline = line.split('|')
            if len(sline) != 
            
        
