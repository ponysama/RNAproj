import os,sys,getopt,re,random
from copy import deepcopy

import dendropy

###############################################################################

def findMatch(struct):
    stk = []
    bp = {}
    for x in range(0,len(struct)):
        if (struct[x] == '('):
            stk.append(x)
        elif (struct[x] == ')'):
            bp[stk.pop()] = x
        else:
            pass
    return bp

def similarityBP(seqDict,struct):
    bp = findMatch(struct)
    rlen = len(seqDict)
    li = list(seqDict.itervalues())
    clen = len(bp)*2
    maxdiff_positions =[]
    tmpdiff = 0
    maxdiff =0
    diff = 0
    for l,r in bp.iteritems():
        for j in range(rlen-1):
            for k in range(j+1,rlen):
                if(li[j][l]!=li[k][l]):
                    diff+=1
                    tmpdiff+=1
                if(li[j][r]!=li[j][r]):
                    diff+=1
                    tmpdiff+=1
        if(tmpdiff>=maxdiff):
            if(tmpdiff>maxdiff):
                maxdiff = tmpdiff
                maxdiff_positions = [(l,r)]
            else:
                maxdiff_positions.append((l,r))
        tmpdiff = 0
    sim = 1 - diff/(float)((((rlen-1)*rlen)/2)*clen)
    sim_difpos = 1- maxdiff/(float)((rlen-1)*rlen) if rlen!=1 else 1
    return [sim,maxdiff_positions,sim_difpos]

def trimLeaves(root,alignsDict,tree):
    for nd in root.postorder_iter():
        if(nd.is_leaf()):
            if (nd.val == {}):
                parentNode = nd.parent_node
                parentNode.remove_child(nd)
            else:
                pass
        else:
            pass
        
        
def trimBranches(root):
    for nd in root.postorder_iter():
        if (len(nd.child_nodes()) == 1):
            if (nd.level() == 0):
                child = (nd.child_nodes())[0]
                grandChildren = [child.remove_child(ch) for ch in child.child_nodes() ]
                nd.remove_child(child)
                for gch in grandChildren:
                    nd.add_child(gch)
            else:
                parentNode = nd.parent_node
                child = (nd.child_nodes())[0]
                nd.remove_child(child)
                parentNode.remove_child(nd)
                parentNode.add_child(child)
        else:
            pass
        
        
# count different loop types

def usage(softname):
    print "%s <arguments>" % softname
    print "-c   Minimum number of clusters"
    print "-C   Maximum number of clusters"
    print "-f   Input file (Required)"
    print "-h   This message"
    print "-m   Minimum MSA size (columns)"
    print "-M   Maximum MSA size (columns)"
    print "-o   Print data in output files (One per cluster of family)"
    print "-s   Maximum number of RNA families (select families randomly if number exceed threshold)"
    print "-S   Maximum number of sequences in RNA families (select sequences randomly if number exceed threshold)"
    print "-x   Select Rfam family (use Rfam ID)"
    sys.exit(1)

###############################################################################

# count different loop types

def ssastats(structure):
    
    stack=[];
    prev_openk  = -1;
    prev_closek = -1;
    sse_counter = {'hairpin':0, 'stack':0, 'bulge':0, 'internal': 0, 'loop': {}, 'bp':0};
    stemstack = []
    
    lindex = range(len(structure));
    for k in lindex:
        carac = structure[k];
        if carac=='(':

            sse_counter['bp']+=1
            stack.append(k);
            if len(stemstack)>0 and prev_closek>=0:
                stemstack[-1]+=1
            stemstack.append(0)
            prev_openk  = -1;
            prev_closek = -1;
        
        elif carac==')':
            closek=k;
            openk=stack.pop()
            nstem = stemstack.pop()
            
            if prev_openk<0 and prev_closek<0:
                sse_counter['hairpin']+=1
            else:
                leftgap  = prev_openk-openk-1
                rightgap = closek-prev_closek-1
                
                if leftgap==0 and rightgap==0:
                    sse_counter['stack']+=1
                else:
                    if nstem==0:
                        if leftgap==0 or rightgap==0:
                            sse_counter['bulge']+=1
                        else:
                            sse_counter['internal']+=1
                    if not sse_counter['loop'].has_key(nstem):
                        sse_counter['loop'][nstem]=0
                    sse_counter['loop'][nstem]+=1
            prev_openk = openk
            prev_closek = closek

    return sse_counter


###############################################################################

def fitstruct2seq(seq,struct,removeNonWC=True):
    
    wcbp = {('G','C'),('C','G'),('A','U'),('U','A'),('G','U'),('U','G')}
    nonwcbp = 0
    
    if len(seq)!=len(struct):
        print "sequence and structure lengths do not match";
        sys.exit(1);
    
    # build bpdic consistent with seq
    bpdic={}
    stack=[]
    for k in range(len(struct)):
        carac=struct[k]
        if carac=='(':
            stack.append(k)
        elif carac==')':
            closei=stack.pop();
            openi=k
            if not (seq[openi],seq[closei]) in wcbp:
                nonwcbp += 1
            if seq[openi]!='.' and seq[closei]!='.':
                bpdic[openi]=closei
                bpdic[closei]=openi
    # build structure
    newseq=''
    newstruct=''
    for k in range(len(seq)):
        nt=seq[k]
        if nt!='.':
            newseq+=nt
            if bpdic.has_key(k):
                if bpdic[k]>k:
                    newstruct+='('
                else:
                    newstruct+=')'
            else:
                newstruct+='.'
    
    return newseq,newstruct,nonwcbp

##########################################################################################

def hashcode(sequenceWithGaps):
    scode=''
    for i in range(len(sequenceWithGaps)):
        if sequenceWithGaps[i]=='.':
            scode += '1'
        else:
            scode += '0'
    return scode

def clusterUnGapped(data):
    clusters={}
    entrymap = {}
    for myid,myseq in data['seqs'].iteritems():
        mycode = hashcode(myseq)
        if not entrymap.has_key(mycode):
            mykey = len(entrymap)
            entrymap[mycode] = mykey
            clusters[mykey] = {}
        else:
            mykey = entrymap[mycode]
        clusters[mykey][myid]=myseq
    return clusters

##########################################################################################

def readRfamDB(filename):
    
    seqline_re = re.compile("(\S+)(\s+)(\S+)");
    
    fh = open(filename,"r")

    info={};
    id=''
    
    for myrow in fh:
        cline = myrow.strip()
        if len(cline)>0:
            if cline[0]=='#':
                if cline.startswith("#=GC SS_cons"):
                    if not info[id].has_key('consensus'):
                        info[id]['consensus']='';
                    info[id]['consensus'] += cline[13:].lstrip().replace('<','(').replace('>',')')
                if cline.startswith("#=GF AC"):
                    id = cline[8:].replace(' ','');
                    info[id]={};
                    info[id]['seqs']={}
                if cline.startswith("#=GF ID"):
                    info[id]['name']=cline[8:].replace(' ','');
                if cline.startswith("#=GF SQ"):
                    info[id]['nseq']=int(cline[8:].replace(' ',''));
            elif cline.startswith("//"):
                if not info[id].has_key('name'):
                    print "WARNING: %s missing name." % (id);
                if not info[id].has_key('nseq'):
                    print "WARNING: %s missing number of sequences." % (id);
                if not info[id].has_key('consensus') or len(info[id]['consensus'])==0:
                    print "WARNING: %s missing consensus structure" % (id);
            else:
                # store sequences
                fields = seqline_re.match(cline);
                if fields:
                    seqid = fields.group(1);
                    seqcontent = fields.group(3);
                else:
                    print "Cannot parse \"%s\"" % cline;
                    sys.exit(1);
                if not info[id]['seqs'].has_key(seqid):
                    info[id]['seqs'][seqid] = '';
                info[id]['seqs'][seqid] += seqcontent;

    fh.close()
                    
    return info

##########################################################################################

def select(ssa,params):

    if len(ssa)<params['min'] or len(ssa)>params['max']:
        return False

    if ssa.count('a') or ssa.count('A'):
        return False

    return True


##########################################################################################

def cleancluster(data,consensus):
    normalizeddata = {}
    for myclusterid,myclusterdata in data.iteritems():
        for seqid,seqdata in myclusterdata.iteritems():
            cleanseq,cleanssa,nonwcbp = fitstruct2seq(seqdata,consensus)
            #kengdie
            #nonwcbp = 0
            mystats = ssastats(cleanssa)
            bpthreshold = len(cleanssa)/5
            # insert sequence with only WC base pairs & at least 20% of base pairs
            if nonwcbp == 0 and mystats['bp']>bpthreshold:
                if not normalizeddata.has_key(myclusterid):
                    normalizeddata[myclusterid] = {}
                    normalizeddata[myclusterid]['sequence'] = {}
                normalizeddata[myclusterid]['sequence'][seqid] = cleanseq
                if not normalizeddata[myclusterid].has_key('consensus'):
                    normalizeddata[myclusterid]['consensus'] = cleanssa
    return normalizeddata


##########################################################################################

def writedataset(data,params):
    
    outfile = data['id'] + '.fasta'
    fh = open(outfile,'w')
    clusteridx = data['largestcluster']
    
    if params['famsize'] > 0 and params['famsize'] < len(data['cluster'][clusteridx]['sequence']):
        idlist = random.sample(list(data['cluster'][clusteridx]['sequence'].keys()),params['famsize'])
    else:
        idlist = list(data['cluster'][clusteridx]['sequence'].keys())
    
    for seqid in idlist:
        seqdata = data['cluster'][clusteridx]['sequence'][seqid]
        fh.write('>'+seqid+'\n')
        fh.write(seqdata+'\n')

    fh.write('>consensus structure\n')
    fh.write(data['cluster'][clusteridx]['consensus']+'\n')
    
    print "%s\t%d/%d" % (data['id'],len(idlist),len(data['cluster'][clusteridx]['sequence']))
                      
    fh.close()

##########################################################################################

def main(params):

    fulldb = readRfamDB(params['infile'])
        
    filterdb = {}
    catalog = {}
    for famid,data in fulldb.iteritems():
        if select(data['consensus'],params):
            # cluster ungapped alignments
            ungappedclusterdata = clusterUnGapped(data)
            # remove sequences with non-canonical bp
            clusterdata = cleancluster(ungappedclusterdata,data['consensus'])
            # find largest cluster
            largestcluster = (-1,-1)
            for clusterid,clusterset in clusterdata.iteritems():
                if len(clusterset['sequence'])>largestcluster[1]:
                    largestcluster = (clusterid,len(clusterset['sequence']))
            #print famid,largestcluster
            # select cluster based on size
            if largestcluster[1] >= params['minclustersize'] and largestcluster[1] <= params['maxclustersize']:
                filterdb[famid] = deepcopy(data)
                filterdb[famid]['id'] = famid
                filterdb[famid]['cluster'] = clusterdata
                filterdb[famid]['largestcluster'] = largestcluster[0]
                mystats = ssastats(filterdb[famid]['cluster'][largestcluster[0]]['consensus'])
                # sort family in according to structure features (i.e. multi-loop)
                maxloopidx=0
                if len(mystats['loop'])>0:
                    for size,count in mystats['loop'].iteritems():
                        if size>maxloopidx:
                            maxloopidx = size
                if not catalog.has_key(maxloopidx):
                    catalog[maxloopidx]=[]
                catalog[maxloopidx].append(famid)

    for famidx,famlist in catalog.iteritems():
        # no family selected
        if params['id']=='':
            if len(famlist)>params['select']:
                fam2write = random.sample(famlist,params['select'])
            else:
                fam2write = famlist
            print '* Category %d: %d/%d' % (famidx,len(fam2write),len(famlist))
            for myfam in fam2write:
                if params['output']:
                    writedataset(filterdb[myfam],params)
                else:
                    #print myfam
                    #kengdie
                    """tmpdata =filterdb[myfam]
                    myoridict = {}
                    clusteridx = tmpdata['largestcluster']
                    idlist = list(tmpdata['cluster'][clusteridx]['sequence'].keys())
                    for seqid in idlist:
                        seqdata = tmpdata['cluster'][clusteridx]['sequence'][seqid]
                        myoridict[seqid] = seqdata
                    simBP = similarityBP(myoridict,tmpdata['cluster'][clusteridx]['consensus'])
                    if(simBP[2]>0.7 and simBP[2]<=0.8):
                        print myfam
                        print "general similarity at all BP sites:" + str(simBP[0])
                        print "least similar BP sites:" + str(simBP[1])
                        print "similarity at the least similar:" +str(simBP[2]) """
                    treefile = "C:\\Users\\zy\\workspace\\rfamPro\\rfam1\\trees\\" + myfam + ".seed_tree"
                    treefileString = (open(treefile).read()).replace('/','#')
                    treefileString = treefileString.replace('\'','')
                    treefileString = treefileString.replace('=','')
                    tree = dendropy.Tree.get_from_string(treefileString, schema="newick")
                    root = tree.seed_node
                    ##############################
                    data = filterdb[myfam]
                    mydict = {}
                    myoridict = {}
                    clusteridx = data['largestcluster']
                    idlist = list(data['cluster'][clusteridx]['sequence'].keys())
                    for seqid in idlist:
                        seqdata = data['cluster'][clusteridx]['sequence'][seqid]
                        mydict[seqid] = seqdata
                        myoridict[seqid] = seqdata
                        #kengdie
                        parts = seqid.split('/')
                        parts2 = parts[1].split('-')
                        seqid2 = parts[0] + '/' + parts2[1] + '-' + parts2[0]
                        mydict[seqid2] = seqdata
                    if (root.distance_from_tip()>=2):
                        
                        for leaf in root.leaf_nodes():
                            seqName = (((leaf.get_node_str()).split())[1]).replace('#','/')
                            
                            if(mydict.has_key(seqName)):
                                leaf.val={'A':0}
                            else:
                                parentNode = leaf.parent_node
                                parentNode.remove_child(leaf)
                        
                        trimLeaves(root,mydict,tree)
                        trimBranches(root)
                        if (root.distance_from_tip()>=3):
                            #print(tree.as_ascii_plot())
                            simBP = similarityBP(myoridict,data['cluster'][clusteridx]['consensus'])
                            if(simBP[2]<0.85):
                                print myfam
                                print simBP[2]
                    
                    
        else:
            if(params['id'] == 'RF01749'):
                    HAYY= 1+1;
            if params['output']:
                if params['id'] in famlist:
                    writedataset(filterdb[params['id']],params)
                else:
                    print >>sys.stderr,"Cannot find family %s" % (params['id'])
            else:
                print params['id']

##########################################################################################

if __name__ == '__main__':
    
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:om:M:c:C:s:S:x:", ["help", "file=", "out", "min=","max=","minC","maxC","select","famsize","id="]);
    except getopt.GetoptError:
        usage(sys.argv[0]);
        
        
    params = {'infile':'C:\Users\zy\workspace\TestProject\rfam\Rfam.seed','output':False,'min':0,'max':1000000,'pk':False,'minclustersize':0,'maxclustersize':1000000,'select':1000000,'famsize':0,'id':''}
        
    argStart=len(sys.argv);
    for o,a in opts:
        if o in ("-h", "--help"):
            usage(sys.argv[0]);
        if o in ("-f", "--file"):
            params['infile'] = a;
            argStart-=2;
        if o in ("-o", "--out"):
            params['output'] = True;
            argStart-=1;
        if o in ("-m", "--min"):
            params['min'] = int(a)
            argStart-=2;
        if o in ("-M", "--max"):
            params['max'] = int(a)
            argStart-=2;
        if o in ("-c", "--minC"):
            params['minclustersize'] = int(a)
            argStart-=2;
        if o in ("-C", "--maxC"):
            params['maxclustersize'] = int(a)
            argStart-=2;
        if o in ("-s", "--select"):
            params['select'] = int(a)
            argStart-=2;
        if o in ("-S", "--famsize"):
            params['famsize'] = int(a)
            argStart-=2;
        if o in ("-x", "--id"):
            params['id'] = a
            argStart-=2;

    main(params);



