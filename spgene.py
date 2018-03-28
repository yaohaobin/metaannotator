import sys
splabel = open(sys.argv[1],'r')

ctggene = open(sys,argv[2],'r')
ctgsp = {}
sphit = []
for line in splabel:
    info = line.split('\t')
	ctgsp[info[0]] = info[1]
	
splabel.close()

for line in genefile:
    info = line.split('\t')
    frag = info[0].split('_')

    gene = info[1]
    ctgid = frag[0]+'_'+frag[1]
    if sp.has_key(ctgid):
        identity = float(info[2])
    if identity >= 70:
        ctgsp = sp[ctgid]
		sphit.append[(ctgsp,gene)]
        if not spgene.has_key(ctgsp):
            spgene[ctgsp] = 1
        else:
            spgene[ctgsp] = spgene[ctgsp] + 1

genefile.close()

spgene = {}.fromkeys(spgene).keys()
critical = {}

for k,v in spgene.items():
    if v>= 10:
	    critical[k] = v
		
print len(critical)
criticalhit = {}
criticalhit ['0'] = []
		
for (sp,gene) in sphit:
    if critical[k].has_key(sp):
        if criticalhit.has_key(sp) and len(criticalhit[sp] <= 20):
		    criticalhit[sp].append(gene)
		else:
		    criticalhit[sp] = [gene]
			
reportlist = []
for k,v in criticalhit.items():
    reportlist.append((k,v))
	
reportlist = sorted(reportlist,reverse=True,key=lambda d:len(d[1]) )

reportlen = len(reportlist)
if reportlen > 10:
    reportlen = 10
for i in range(0,reportlen):
    print reportlist[i][0]
	for gene in reportlist[i][1]:
	    print gene,
	print ''

	    


