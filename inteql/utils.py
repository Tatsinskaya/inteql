

# Get resolution bin of a region
def position2matrixBin(s): return int(s) - (int(s) % 5000)

# Rewriting features
def variantId2chr(s): return 'chr'+ s.split('_')[0]
def variantId2chrNum(s): return s.split('_')[0]
def variantId2pos(s): return s.split('_')[1]
def variantId2end(s): return int(s.split('_')[1]) + len(s.split('_')[3])
def geneIdVersion2geneId(s): return s.split('.')[0]

def regionId2pos(s): return int(s.split('|')[1].split(':')[1].split('-')[0])

def addChrPrefix(s): return 'chr' + str(s)

def split_list(s): return s.split(',')