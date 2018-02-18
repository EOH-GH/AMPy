import os
import sys
if sys.version_info[0]>=3:
    import subprocess
else: 
    import commands

import re ## for substitutions
from collections import Counter
import datetime
import glob
import pandas
import json
import numpy
import math

#folderTped='~/AMPy'
#tpeds=glob.glob(folderTped+'/*tped')
#tfams=[re.sub('.tped','.tfam',i) for i in tpeds]  
#folderFiles=input("Enter path to references: ")
#tped=input("Enter path to tped: ")
#All00pth=input("Enter path to All00: ")
#folderFiles=input("Enter path of folder containing Arithmetic Medoid references: ")

##Paths to tped/references. Set wd.
os.chdir(os.getcwd())
tped=input("Enter full path to tped: ")
folderFiles=[input("Enter path to references: ")]
#folderFiles=['/home/eoghan/2015-Paper-Thesis']
#folderFiles=['/home/eoghan/SonicMedoids','/home/eoghan/SonicMedoids','/home/eoghan/SonicMedoids/1k','/home/eoghan/SonicMedoids/chipnone','/home/eoghan/SonicMedoids/hell']
#tped='~/AMPy/BrPolAdam.tped'


#files=glob.glob(folderFiles+'/medoidArithmetic*ods')

##Assemble list of reference paths
files=[]
for z in folderFiles: files=files+glob.glob(z+'/medoidArithmetic*ods')

print(files)

##Read references, import pop names, rsids and genotypes into lists
popnams=[]
rsidls=[]
genols=[]
for i in files:
    #print(str(i))
    nam=i.split("_")[1]+'.'+i.split("_")[2]
    popnams.append(nam)
    f=open(i)
    lines=f.readlines()
    f.close()
    lines[1] = lines[1].split(" ")
    del lines[1][0]
    lines[0] = lines[0].split(" ")
    del lines[0][0]
    #print(len(lines[0]))
    if i==files[0]:
        print('\033[1m'+'Loading arithmetic references.'+'\033[0m')
        sharedsnps=lines[0]
        sharedsnps=set(sharedsnps)
    if i!=files[0]: sharedsnps=sharedsnps & set(lines[0])
    #print('Reading '+nam+': '+str(len(sharedsnps))+' shared SNPs between '+nam+' and previous loaded references.')
    print('Reading '+'\033[1m'+nam+'\033[0m'+': '+'\033[1m'+str(len(sharedsnps))+'\033[0m'+' shared SNPs between '+nam+' and previous loaded references.')
    #if i!=files[0]: sharedsnps=sharedsnps & set(lines[0])
    lines[1]=[re.sub('\n$','',i) for i in lines[1]]
    ##cut off last 3 characters for all entries
    for i in range(0,len(lines[0])): lines[0][i]=lines[0][i][:-3]
    lines[0][len(lines[0])-1]=lines[0][len(lines[0])-1][:-1]
    rsidls.append(lines[0])
    lines[1]=list(map(float,lines[1]))
    genols.append(lines[1])


sharedsnps2=list(sharedsnps)

##Take dbsnp orientation from reference files
dbstrd={}
for i in range(0,len(sharedsnps2)-1):
        #dbstrd.append(sharedsnps2[i][-2:])
        dbstrd[sharedsnps2[i][:-3]]=sharedsnps2[i][-2:]
        sharedsnps2[i]=sharedsnps2[i][:-3]


##SNPs shared between references
print('Number of shared SNPs between references used: '+'\033[1m'+str(len(sharedsnps))+'\033[0m')

##Read tped, remove all columns except genotypes
df=pandas.read_csv(tped,header=None,sep=' ')
del df[0];del df[2];del df[3]#;del df[1]
df.index=df[1]
dfsnps=df[1]
del df[1]


##Check that there are common SNPs between reference and sample
if len(list(set(dfsnps) & set(sharedsnps2)))==0:
    print('Error, no shared SNPs between sample and all References used. References use the SNP ids from dbSNP, are the sample SNPs annotated the same way? Is the file submitted of the correct (Tped) format? Another possible solution is to limit references to those using a similar SNP panel or made from whole genome sources in the case of data not from SNP panels.')
    exit()

df=df.loc[set(dfsnps) & set(sharedsnps2)]
dfsnps=list(df.index)

##Read tfam
dfid=pandas.read_csv(re.sub('.tped','.tfam',tped),header=None,sep=' ')
dfid=dfid[1]

##Scan tped for all unique characters for each loci, put them in a list, detect any with more than 3 alpha chars as these are indels or errors.
uniloc=[]
indelrs=[]
for i in list(df.index):
       uniloc.append(''.join(set(df.loc[i][1:])))
       indels=map(lambda x:len(x)>2 and x.isalpha()==True,df.loc[i][1:])
       #if any True in indels: indelrs.append(i)
       if len(set(indels))==2: indelrs.append(i)

       
##Print warning on indels.
if len(indelrs)>0:
    if lst.count==1:print(str(len(indelrs))+' SNP with 2 or more character entries detected, removing, these could be indels or sequencing/format errors.')
    if lst.count>1:print(str(len(indelrs))+' SNPs with 2 or more character entries detected, removing, these could be indels or sequencing/format errors.')


#df2=df.iloc[:,range(4,len(df.columns)+3,2)]
#example=['CG','TA','TTA','0GT']
#map(lambda x:True if x!='CG' or x!='AT' or x.isalpha()==True or len(x)==2 else False,example)

#!(map(lambda x:True if ''x=='CG' or x=='AT' or '0' in x or len(x)>2 else False,example))
#examplers=['rs1','rs2','rs3','rs4']
#len(example.isin(['CG']))
       

##Flag warning for missingness.
missmap=set(list(map(lambda x:'0' in x,uniloc))) if sys.version_info[0]>=3 else set(map(lambda x:'0' in x,uniloc))

if len(missmap)>1:
       #lst=map(lambda x:'0' in x,uniloc)
       lst=list(map(lambda x:'0' in x,uniloc)) if sys.version_info[0]>=3 else map(lambda x:'0' in x,uniloc)
       if lst.count(True)==1:print(str(lst.count(True))+' SNP with Missingness detected, removing. If minor in number and loss is unacceptable, try imputation on individuals or removal of individuals with missingness.')
       if lst.count(True)>1:print(str(lst.count(True))+' SNPs with Missingness detected, removing. If minor in number and loss is unacceptable, try imputation on individuals or removal of individuals with missingness.')

##Flag warning of CGs and ATs.
if 'CG' or 'AT' in uniloc:
       lst=list(map(lambda x:x=='CG' or x=='AT',uniloc)) if sys.version_info[0]>=3 else map(lambda x:x=='CG' or x=='AT',uniloc)
       if lst.count(True)==1:print(str(lst.count(True))+' CG/AT Present, strand reading issue, removing.')
       if lst.count(True)>1:print(str(lst.count(True))+' CGs/ATs Present, strand reading issue, removing.')


##Flag any indels
indlmap=set(list(map(lambda x:len(x)>2,uniloc))) if sys.version_info[0]>=3 else set(map(lambda x:len(x)>2,uniloc))
#if len(set(map(lambda x:len(x)>2,uniloc)))>1:
if len(indlmap)>1:
       lst=list(map(lambda x:len(x)>2,uniloc)) if sys.version_info[0]>=3 else map(lambda x:len(x)>2,uniloc)
       if lst.count(True)==1:print(str(lst.count(True))+' SNP with 2 or more character entries detected, removing, these could be indels or sequencing/format errors.')
       if lst.count(True)>1:print(str(lst.count(True))+' SNPs with 2 or more character entries detected, removing, these could be indels or sequencing/format errors.')
       


##Keep only SNPs that are not CG, AT, have no missingness and have no indels.
#rwkp=list(map(lambda x:True if x!='CG' or x!='AT' or x.isalpha()==True or len(x)<=2 else False,uniloc)) if sys.version_info[0]>=3 else map(lambda x:True if x!='CG' or x!='AT' or x.isalpha()==True or len(x)<=2 else False,uniloc)
rwkp=list(map(lambda x:True if x!='CG' and x!='AT' and x.isalpha()==True and len(x)<=2 else False,uniloc)) if sys.version_info[0]>=3 else map(lambda x:True if x!='CG' and x!='AT' and x.isalpha()==True and len(x)<=2 else False,uniloc)
#rwkp.count(True)
df=df.loc[rwkp]


       
##Dealing with Mono CGs
#monoCGs=list(map(lambda x:True if x=='C' or x=='G' else False,uniloc)) if sys.version_info[0]>=3 else map(lambda x:True if x=='C' or x=='G' else False,uniloc) 
#monoCGs=list(map(lambda x:uniloc.index(x) if x=='C' or x=='G' else False,uniloc)) if sys.version_info[0]>=3 else map(lambda x:uniloc.index(x) if x=='C' or x=='G' else False,uniloc) 
monoCGs=[]
for i in range(0,len(dfsnps)):
    if uniloc[i]=='G' or uniloc[i]=='C': monoCGs.append(dfsnps[i])

        
if len(list(set(df.index) & set(sharedsnps2)))==0:
            print('Error, no shared SNPs between sample and all References used. References use the SNP IDs from dbSNP, are the sample SNPs annotated the same way? Another possible solution is to limit references to those using a similar SNP panel or made from whole genome sources in the case of data not from SNP panels.')


##Replace values for MonoCGs
if len(monoCGs)>0:
    for i in monoCGs:
        df.loc[i]=df.loc[i].replace({min(dbstrd[i]): 1}, regex=True)
        df.loc[i]=df.loc[i].replace({max(dbstrd[i]): 2}, regex=True)


##Replace all the remaining As and Ts
df=df.replace({'A': 1}, regex=True)
df=df.replace({'T': 2}, regex=True)

##Remove str() later, likely a more efficient way to do this.
for i in df.index:
    #if 'C' in str(df.loc[i]) and min(df.loc[i])==1: df.loc[i]=df.loc[i].replace({'C': 2},regex=True)
    #if 'C' in str(df.loc[i]) and min(df.loc[i])==2: df.loc[i]=df.loc[i].replace({'C': 1},regex=True)
    #if 'G' in str(df.loc[i]) and min(df.loc[i])==1: df.loc[i]=df.loc[i].replace({'G': 2},regex=True)
    #if 'G' in str(df.loc[i]) and min(df.loc[i])==2: df.loc[i]=df.loc[i].replace({'G': 1},regex=True)
    if 'C' in list(df.loc[i]) and 1 in list(df.loc[i]): df.loc[i]=df.loc[i].replace({'C': 2},regex=True)
    if 'C' in list(df.loc[i]) and 2 in list(df.loc[i]): df.loc[i]=df.loc[i].replace({'C': 1},regex=True)
    if 'G' in list(df.loc[i]) and 1 in list(df.loc[i]): df.loc[i]=df.loc[i].replace({'G': 2},regex=True)
    if 'G' in list(df.loc[i]) and 2 in list(df.loc[i]): df.loc[i]=df.loc[i].replace({'G': 1},regex=True)


##Add both allele scores for each individual.
for i in list(df.index):
       #fst=df.loc[i][range(4,len(df.columns)+3,2)]
       fst=df.loc[i][list(range(4,len(df.columns)+3,2))] if sys.version_info[0]>=3 else df.loc[i][range(4,len(df.columns)+3,2)]
       #snd=df.loc[i][range(5,len(df.columns)+4,2)]
       snd=df.loc[i][list(range(5,len(df.columns)+4,2))] if sys.version_info[0]>=3 else df.loc[i][range(5,len(df.columns)+4,2)]
       genscr=[x+y for x, y in zip(fst, snd)]
       df.loc[i,range(4,len(df.columns)+3,2)]=genscr


df=df[list(range(4,len(df.columns)+3,2))] if sys.version_info[0]>=3 else df[range(4,len(df.columns)+3,2)]
df.columns=dfid

matchset=lambda a, b: [b.index(x) for x in a]
#matchset(set(dfsnps) & set(sharedsnps2),sharedsnps2)
#matchset(set(dfsnps) & set(sharedsnps2),list(dfsnps))

#df.loc[set(dfsnps) & set(sharedsnps2)]
#indsl=df.loc[list(set(dfsnps) & set(sharedsnps2)), 'P55_QIA_CONC']
#popgeno=list(map(float,popgeno))
#snpscore=math.sqrt(sum([(x-y)*(x-y) for x, y in zip(list(indsl), popgeno)]))/len(popgeno)

##Calculate and collect all C values for each individual against each ref
Cvalnams=[]
Ivalnams=[]
for i in range(0,len(popnams)):
    snpscore=[]
    ##
    Cvalnams.append('C_'+popnams[i])
    Ivalnams.append('I_'+popnams[i])
    popgeno=[genols[i][z] for z in matchset(set(dfsnps) & set(sharedsnps2),sharedsnps2)]
    #df.loc[matchset(set(sharedsnps2) & set(dfsnps),list(dfsnps)):, 'P55_QIA_CONC']-popgeno
    #df.loc[list(set(sharedsnps2) & set(dfsnps)):, 'P55_QIA_CONC']-popgeno
    for z in list(dfid):
        #indsl=df.loc[list(set(dfsnps) & set(sharedsnps3)), z]
        indsl=df.loc[:, z]
        snpscore.append(math.sqrt(sum([(x-y)*(x-y) for x, y in zip(list(indsl), popgeno)]))/len(popgeno))
    
    if i==0: AMid=pandas.DataFrame(snpscore)
    if i!=0: AMid=pandas.concat([AMid, pandas.DataFrame(snpscore)], axis=1)



AMid.columns=Cvalnams
AMid.index=dfid
AMid.index.names=['UNIQID']

##Calculate normalised values.
normal100=lambda a: [100-((x-minvec)/(maxvec-minvec)*100) for x in a]
for i in range(0,len(dfid)):
    #colsli=list(AMid.iloc[:, i])
    colsli=list(AMid.iloc[i])
    minvec=min(colsli)
    maxvec=max(colsli)
    normvals=normal100(colsli)
    #if i!=0: normvals=normvals[0:len(normvals)/2]
    if i==0: AMid2=pandas.concat([AMid, pandas.DataFrame(data=normvals,index=Ivalnams,columns=[dfid[i]]).T], axis=1).reindex(AMid.index)
    if i!=0: AMid2.update(pandas.DataFrame(data=normvals,index=Ivalnams,columns=[dfid[i]]).T)

    
print('Number of SNPs shared between '+tped.split("/")[-1]+' and all references after QC: '+str(len(set(dfsnps) & set(sharedsnps2))))

print('Results are output to: '+re.sub('.tped','',tped.split("/")[-1])+'_Inds-'+str(len(dfid))+'_'+'SNPs-'+str(len(set(dfsnps) & set(sharedsnps2)))+'.amid')

##Output results to .amid file.
AMid2.to_csv(re.sub('.tped','',tped.split("/")[-1])+'_Inds-'+str(len(dfid))+'_'+'SNPs-'+str(len(set(dfsnps) & set(sharedsnps2)))+'.amidPy2', sep='\t')

