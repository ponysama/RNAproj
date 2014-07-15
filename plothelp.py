"""inpath = "E:\\M1\\PROJ\\trying\\point7to8.txt"
infile = open(inpath)
outpath = "E:\\M1\\PROJ\\trying\\point7to8_stripped.txt"
outfile = open(outpath,"w")
for line in infile:
    if (line[:2] == "RF"):
        outfile.write(line)
        
outfile.close()

infile.close()
"""
import subprocess
import sys
import dendropy
import re
import matplotlib.pyplot as plt
import numpy as np
import scipy.misc as misc 
from matplotlib.backends.backend_pdf import PdfPages
pp = PdfPages('E:\\M1\\PROJ\\trying\\classical_dot.pdf')

#sys.path = "C:\\Users\\zy\\Downloads\\DendroPy-3.12.0\\dendropy"

"""inpath = "E:\\M1\\PROJ\\trying\\4gengraph\\names.txt"
ifile = open(inpath)
prefix = "-f C:\\Users\\zy\\workspace\\TestProject\\rfam\\Rfam.seed"
for line in ifile:
    famid = line[:7]
    opt = " -x " + famid
    with open("E:\\M1\\PROJ\\trying\\classical\\fa_files\\level4_"+famid+".txt",'w') as outfile:
        subprocess.call(["python","C:\\Users\\zy\\workspace\\rfamPro\\rfam1\\6_10.py","-fC:\\Users\\zy\\workspace\\TestProject\\rfam\\Rfam.seed","-x",famid],stdout=outfile)"""

#with open("E:\\M1\\PROJ\\trying\\lepoint6\\result_"+"RF00666_propagate"+".txt",'w') as outfile:
#    subprocess.call(["python","C:\\Users\\zy\\workspace\\rfamPro\\rfam1\\6_10.py","-fC:\\Users\\zy\\workspace\\TestProject\\rfam\\Rfam.seed","-x","RF00666"],stdout=outfile)
    
def plotEnergy(pathin,fam):
    if fam=='RF01034':
        jj=1
    infile = open(pathin)
    Edict = {}
    lines = infile.readlines()
    for i in range(len(lines)):
        if (lines[i][0]=='>'):
            level = int(lines[i][1])
            
            nl = lines[i+2]
            m = re.search('(?<=-)\w+...', nl)
            strE ='-' + m.group(0)
            numE = float(strE)
            if Edict.has_key(level):
                (Edict[level]).append(numE)
            else:
                Edict[level]=[numE]
    print Edict
    avgDict = {}
    for key,val in Edict.iteritems():
        avgDict[key]= sum(val)/(float)(len(val))
    print avgDict
    li1 = [i for i in avgDict.iterkeys()]
    vals = [avgDict[i] for i in li1]
    
    plt.plot(li1,vals)
    plt.ylabel('energy')
    plt.xlabel('number of generations from root ---'+fam+'(with structure)')
    
    #if fam=='RF01034':
    #    plt.show()
    #axes = plt.subplot(111)
    #xmajorLocator   = plt.MultipleLocator(0.5)
    #ymajorLocator   = plt.MultipleLocator(0.5)
    #axes.xaxis.set_major_locator(xmajorLocator) 
    #axes.yaxis.set_major_locator(ymajorLocator) 
    plt.savefig(pp, format='pdf')
    plt.clf()
#inpath = "E:\\M1\\PROJ\\trying\\RF00666_propagate.txt"
#plotEnergy(inpath)

inpath = "E:\\M1\\PROJ\\trying\\4gengraph\\names.txt"
names = open(inpath)
lines =names.readlines()
for i in range(len(lines)):
    fam = lines[i][:7]
    tmppath = "E:\\M1\\PROJ\\trying\\classical\\med\\lev4" + fam+".txt"
    plotEnergy(tmppath,fam)
pp.close()


"""li = [0,1,2]
li2= [-38.599999999999994,-38.599999999999994,-40.15]
plt.plot(li,li2)
plt.clf()
li3=[0,1,2]
li4=[-50,-80,-90]
plt.plot(li3,li4)
plt.show()"""
