#Zachary Hinz, Abby Smith, Riki Imai, Maya Singapuri
import random
from collections import defaultdict

class ibm2:

 def __init__(self, e_sent, f_sent, iterations, threshhold, alignmentsflag):
     engFile, forFile= open(e_sent, "r"), open(f_sent, "r")


     # P(f|e), alignment_dict
     # P(e, f)/count(e, f), count (e)
     pfe, alignment = {}, {}
     engSent, forSent, fileLength = [], [], 0
  
     # read english File into arraylist
     for line in engFile:
         engSent.append("NULL "+line)
         #engSent.append(line)
         fileLength += 1


     # read foreign file into arraylist
     for line in forFile:
         forSent.append(line)

     engFile.close()
     forFile.close()

     engVocab = set() #added to make IBM model 2 easier
     forVocab = set()

     # set up pfe tables, assigning them a value of 0.01
     for i in range(fileLength):
         for e in list(set(engSent[i].split())):
             engVocab.add(e)
             for f in list(set(forSent[i].split())):
                 forVocab.add(f)   
                 if(e not in pfe):
                     pfe[e] = {f:0.01}
                 elif(f not in pfe[e]):
                     pfe[e][f] = 0.01
                
     print("PFE Table Created")

     # complete alignment alg for given iterations
     for times in range(iterations):
         # instantiate the count tables
         ce = {}
         cef = {}
         for i in range(fileLength):
             # iterate through each pairing
             if(i%1000==0):
                 print("Line Num: "+str(i))
             engWords = list(set(engSent[i].split()))
             for e in engWords:
                 # add word to ce and cef tables as needed
                 if e not in ce:
                     ce[e] = 0
                     cef[e] = {}
                 forWords = list(set(forSent[i].split()))
                 cache = {}
                 for f in forWords:
                     # add word to cef[e] as needed
                     if f not in cef[e]:
                         cef[e][f] = 0
                     # cache to reduce calculations
                     if f in cache:
                         denom = cache[f]
                     else:
                         denom = 0
                         for e_s in engWords:
                             denom+=pfe[e_s][f]
                         cache[f] = denom
                  
                     pFtoE = pfe[e][f]/denom
                     cef[e][f] += pFtoE
                     ce[e] += pFtoE 


         # update p(f->e) after setup        
         for eng in cef:
             for foreign in cef.get(eng):
                 pfe.get(eng)[foreign] = (cef.get(eng)[foreign])/(ce[eng])

         print("Iteration: "+str(times))

     # ***** IBM MODEL 2 *****

     print("start IMB 2")

     # line 2 of pseudocode -- initalize a(i|j, le, lf) to 1/(lf+1) for all i, j, le, lf, goes [f][e][le][lf]
     alignModel2 = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: float))))
     for i in range(0, len(engSent)):
         curEngSent = engSent[i].split()
         curForSent = forSent[i].split()
         lf = len(curForSent)
         le = len(curEngSent)
         value = 1 / (lf)

         for j in range(0, lf):
             for k in range(0, le):
                 f = curForSent[j]
                 e = curEngSent[k]
                 alignModel2[f][e][le][lf] = value

     for times in range(iterations):

        #1 things, line 6 pseudo
        #2 things, line 5 pseudo
        countEF = defaultdict(lambda: defaultdict(float))
        totalF = defaultdict(float)
        # line 9 of pseudocode
        totalE = defaultdict(float)

        #4 things, line 7
        #3 things, line 8
        countAlign = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0.0))))
        totalAlign = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0.0)))

        '''
        for i in range(0, len(engSent)):
            curEngSent = engSent[i].split()
            curForSent = forSent[i].split()

            lf = len(curForSent)
            le = len(curEngSent)
            for f in range(0, lf):
                for e in range(0, le):
                    j = forSent[i].split()[f]
                    k = engSent[i].split()[e]
                    if(j not in countAlign):
                        countAlign[j] = {}
                    if(k not in countAlign[j]):
                        countAlign[j][k] = {}
                    if(le not in countAlign[j][k]):
                        countAlign[j][k][le] = {}
                    if(lf not in countAlign[j][k][le]):
                        countAlign[j][k][le][lf] = {}
                    countAlign[j][k][le][lf] = 0
                if(j not in totalAlign):
                    totalAlign[j] = {}
                if(le not in totalAlign):
                    totalAlign[j][le] = {}
                if(lf not in totalAlign):
                    totalAlign[j][le][lf] = 0
            


            #initializations
        
        for i in range(fileLength):
            for f in list(forSent[i].split()):
                totalF[f] = 0
                for e in list(engSent[i].split()):
                    if(e not in countEF):
                        countEF[e] = {}
                    countEF[e][f] = 0
        '''

        for i in range(fileLength):
            curEngSent = engSent[i].split()
            curForSent = forSent[i].split()

            le = len(curEngSent)
            lf = len(curForSent)

            # line 11 pseudocode -- compute normalization
            for j in range(0, le):
                e = curEngSent[j]
                totalE[e] = 0
                for k in range(0, lf):
                    f = curForSent[k]
                    temp = pfe[e][f]*alignModel2[f][e][le][lf]
                    totalE[e] += pfe[e][f] * alignModel2[f][e][le][lf] # line 15

            # line 18 pseudocode -- collect counts 
            for j in range(0, le):
                e = curEngSent[j]
                for k in range(0, lf):
                    f = curForSent[k]
                    delta = pfe[e][f] * alignModel2[f][e][le][lf] / totalE[e]
                    countEF[e][f] += delta
                    totalF[f] += delta # line 23
                    countAlign[f][e][le][lf] += delta
                    totalAlign[f][le][lf] += delta

        #line 29 pseudo -- estimate probabilites

        '''
        # t(e|f) = 0 for all e, f -- line 30
        for e in engVocab: # collected "vocab" when reading in the sentences
            for f in forVocab:
                pfe[e][f]= 0
        
        # a(i| j, le, lf) = 0 for all i, j, le, lf -- line 31
        for i in alignModel2:
            for j in alignModel2[i]:
                for le in alignModel2[i][j]:
                    for lf in alignModel2[i][j][le]:
                        alignModel2[i][j][le][lf] = 0
        '''

        pfe[e][f] = defaultdict(lambda: defaultdict(lambda: 0.0))
        alignModel2[i][j][le][lf] = defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: defaultdict(lambda: 0.0))))

        # line 32 & 33, new lexical translation probabilities
        for e in engVocab:
            for f in forVocab:
                #if e in pfe and f in pfe[e] and e in countEF and f in countEF[e] and f in totalF:
                pfe[e][f] = countEF[e][f]/totalF[f]

        # line 35, new alignment probabilities
        for i in range(0, fileLength):
            curEngSent = engSent[i].split()
            curForSent = forSent[i].split()
            lf = len(curForSent)
            le = len(curEngSent)
            for j in range(0, lf):
                for k in range(0, le):
                    f = curForSent[j]
                    e = curEngSent[k]
                    alignModel2[f][e][le][lf] = countAlign[f][e][le][lf] / totalAlign[f][le][lf]



        # for all i, j, le, lf do a(i| j, le, lf) = countA(i|j, le, lf)/totalA(i|j, le, lf)
        

     print("pfe")
     #print(pfe)


     # print p(f|e) table sorted by english and foreign words
     # print("\n"+e_sent+" "+f_sent+" p(f|e) Table")
  
     if fileLength > 10:
         rand10 = {}
         for _ in range(30):
             eng, forChoices = random.choice(list(pfe.items()))
             rand10[eng] = forChoices
         for eng, forWords in sorted(rand10.items(), key=lambda x: x[0]):
             for key, val in sorted(forWords.items(), key=lambda x: x[0]):
                 if val > threshhold:
                     print(eng+"\t"+key+"\t"+str(val))
     else:
         for eng, forWords in sorted(pfe.items(), key=lambda x: x[0]):
             for key, val in sorted(forWords.items(), key=lambda x: x[0]):
                 highestProb = 0
                 if val > threshhold:
                     print(eng+"\t"+key+"\t"+str(val))
      

     # generate alignments if needed
     if(alignmentsflag == True):
         picks = []
         if(fileLength > 10):
             # pick 10 random sentences
             for _ in range(10):
                 picks.append(random.randint(0, fileLength-1))
         else:
             for i in range(fileLength):
                 picks.append(i)
      
         for idx in picks:
             print("\n"+"Foreign Sentence: "+forSent[idx]+"\n"+"English Sentence: "+engSent[idx])
             # for each foreign word
             for foreign in forSent[idx].split():
                 # go through all english words in sentence
                 ForToEngidx = -1
                 highestProb = 0
                 engWords = engSent[idx].split()
                 ForToEngWord = "UNK"
                 for j in range(len(engWords)):
                     # pick max of all the p(f|e) vals
                     if(pfe.get(engWords[j])[foreign] > highestProb):
                         highestProb = pfe.get(engWords[j])[foreign]
                         ForToEngidx = j
                         ForToEngWord = engWords[j]
                 print(foreign+" -> "+ForToEngWord)


# requires 5th arg (boolean) stating if you want alignments for this corpus
#test0 = em('test.en 1', 'test.es 1', 20, 0.1, True)
#test1 = em('test.en', 'test.fr', 10, 0, True)
#fullRun = em('es-en.10k.en', 'es-en.10k.es', 10, 0.3, True)
#test2 = ibm2("test.en", "test.fr", 4, 0.3, True)
fullRun = ibm2('es-en.1k.en', 'es-en.1k.es', 5, 0.3, True)
