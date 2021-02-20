def convert(s): 
        new = "" 
        for x in s: 
            new += x  
        return new



import pandas as pd
import matplotlib.pyplot 
from matplotlib.pyplot import gca

filename='cores' 
file = open(filename,mode='r')
file_read = file.read()
file.close()
s = file_read
 
#converting everything to lower case
if s.islower() == False:
    s=s.lower()

#finding number of sequences in the input file filename
count_seq = s.count('>')
 
#splitting sequences to individual arrays
s= s.split('>')

#separating sequences from names
arr_seq=['']* count_seq
for i in range(1, len(s)):
    
    sm = str(s[i])
    lines=  sm.count('\n')
    for ix in range(1,lines):
        arr_seq[i-1] += sm.split('\n')[ix]

    
arr_seqname=[0]* count_seq
for i in range(1, len(s)):
    sm = str(s[i])
    arr_seqname[i-1] = sm.split('\n')[0]
    
 
name = arr_seqname 
s_array = arr_seq
seqs= []

def find_vdre(s,filen):
    for im in range (0,count_seq): 
        h=open('Output_%s_%s'%(filen,arr_seqname[im]),'w')
        s= s_array[im]
        
        if filen =='rev':
            s = list(s)
            l = len(s) 
            
            for j in range (0,l):
                if (s[j]=='g'):
                    s[j] = 'c'
                elif (s[j]=='a'):
                    s[j] = 't'    
                elif (s[j]=='c'):
                    s[j] = 'g'
                elif (s[j]=='t'):
                    s[j] = 'a'
             
            s2 = s[::-1]    
            s2 = convert(s2 )    
            s = str(s2) 
          
        
        l = len(s)  
        s3 = ""
        m2 =""
        m1 =""
        a1 =""
        a2 = ""
        b1 =""
        b2 =""
        c1=""
        c2=""
        d1=""
        d2=""
        
        count_all=0
        for i in range (0, l-13):
            
            for n in range (8,9):
     
                    count = 0
                    if ((s[i] == s[i+n+1]) and (s[i]== 'g')):
                        if ((s[i-1]== 'g') or (s[i-1]== 'a')): 
                            count+=1        
                        if ((s[i+1]== 'g') or (s[i+1]== 'a')or (s[i+1]== 't')): 
                            count+=1            
                        if ((s[i+2]== 'g') or (s[i+2]== 't')): 
                            count+=1                
                        if ((s[i+3]== 'c') or (s[i+3]== 't')):
                            count+=1                     
                        if ((s[i+4]== 'g') or (s[i+4]== 'a')): 
                            count+=1
                    count2=0
                    if (count > 3):   
                        if  ((s[i+n]== 'g') or (s[i+n]== 'a')): 
                            count2+=1        
                        if ((s[i+n+2]== 'g') or (s[i+n+2]== 'a')or (s[i+n+2]== 't')): 
                            count2+=1            
                        if ((s[i+n+3]== 'g') or (s[i+n+3]== 't')): 
                            count2+=1                
                        if ((s[i+n+4]== 'c') or (s[i+n+4]== 't')):
                            count2+=1                     
                        if ((s[i+n+5]== 'g') or (s[i+n+5]== 'a')): 
                            count2+=1
                    if (count2 > 3):    
                           count_all+=1 
    
                           m1 = i-1
                           m2 = m2 + str(m1) + '--'
                           
                           a1 = i + n +6
                           a2 = a2 + str(a1) + '--'
                           
                           b1 = n-5
                           b2 = b2 + str(b1) + '--'
                           
                           c1 = count +1
                           c2 = c2 + str(c1) + '--'
                           
                           d1 = count2+1
                           d2 = d2 + str(d1) + '--'
                           
                           sp3 = s[i-1:i+n+6]
                           s3 = s3 + sp3  + '--'                  
        
        seq = s3.split('--')
        on = m2.split('--')
        off =a2.split('--')
        seq_type=b2.split('--')
        correct1=c2.split('--')
        correct2 =d2.split('--')
        
        seq.remove('')
        on.remove('')
        off.remove('')
        seq_type.remove('')
        correct1.remove('')
        correct2.remove('')
        
         
        for i in range(0, len(on)):
            c = int (on[i])
            on[i] = c
            c2 = int (off[i])
            off[i] = c2
            d= int (seq_type[i])
            seq_type[i] =d
         
    
        n = 25
        n=int(n)
        
        if (l%n == 0):
            ls = l/n
        else:
            ls = (l/n) + 1
            
        ls = int(ls) * n
        count_graph = [0] * (ls)
        for i in range (0,ls):
            for j in range (0,len(on)):
                 
                    if(i<=on[j]<=i+n)  and (i<=off[j]<=i+n) :
                        count_graph[i] +=1
     
        x_smooth = list(range(0,len(count_graph)))
        
        x_axis = x_smooth
        y_axis = count_graph
        
        image= matplotlib.pyplot.figure(figsize=(3.3464566929,0.826))
        matplotlib.pyplot.xlim(0,l) 
        matplotlib.pyplot.ylim(0,6)
        #matplotlib.pyplot.xlabel('Nucleotide Position on HBV Genome')
        #matplotlib.pyplot.ylabel('Number of VDREs')
        matplotlib.pyplot.xticks([])
        matplotlib.pyplot.yticks([0,2,4,6])
        matplotlib.pyplot.yticks(fontname = "Arial", fontsize ='11') 
        matplotlib.pyplot.bar(x_axis,y_axis, linewidth = 25) 
        image.savefig('_%s_%s'%(filen,arr_seqname[im]), bbox_inches='tight', pad_inches=0, dpi = 400) #'''
        
        
        for i in range(len(on)):
            h. write(str(i)+ '    '+ str(seq [i]) + '    ' + str(on[i])+'    '+str(off[i]) + '\n')
        h.close()
        
        if len(seq) >0:
            df = pd.read_csv('Output_%s_%s'%(filen,arr_seqname[im]), delimiter= '    ')
            df.columns = ['index','seq', 'start', 'stop'] 
            df.to_csv('VDREs_%s'%(arr_seqname[im])+'.csv',index=False)  
    return df      
            
get_vdre=find_vdre(s,'fwd')
get_vdre_rev=find_vdre(s, 'rev')
     
         
           
