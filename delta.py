########## COMPUTING DELTA-STATISTICS FROM PLINK TRAW-FILES #########################################################

import math

# mymean
# Computes half (because diploid) the mean of the provided list.
#
# list	A list of strings that are allele counts of diploid individuals from a population or missing data coded as "NA".
# 
# Returns the allele frequency unless all data was missing, in which case returns "NA".
def mymean(list):
    sum = 0
    amo = 0
    for item in list:
        if item != "NA":
            sum = sum + int(item)
            amo = amo + 1
    if amo == 0:
        return("NA")
    else:
        return(0.5*sum/amo)

# Delta
# Computes estimators of Delta-statistics and variances of those estimators.
#
# traw		A PLINK traw-file created with the command --recode A-transpose. Contains a header row and one row for each genetic locus.
#		First six columns are chromosome, ID, position in cmorgans, base-pair coordinate, counted allele and the reference allele.
#		Of these only the chromosome and the base-pair coordinate is used by the function.
#		The rest of the columns are allele counts of diploid individuals, missing data coded as "NA".
# populations	List of K lists. The K:th contain indices of the columns corresponding to the K:th population.
#		Keep in mind Python starts counting from zero and the first six columns contain other stuff, so the smallest possible index is 6.
# left		List of patterns that contribute positively to the Delta-statistic, written as words in {A, B}^K where B is the counted allele.
# right		List of patterns that contribute negatively to the Delta-statistic, written as words in {A, B}^K where B is the counted allele.
# blocksize	The size of blocks used, defaults at two hundred thousands (if the parameter correction is set to True, increase this to at least two millions).
#		A new block starts when the chromosome changes or the blocksize is reached in base pair coordinates.
#               To suppress the use of blocks (only when the parameter correction is False), set the blocksize to 0.
# correction	Parameter that decidces what to do with the designated blocks.
#               The default option is False: window-wise computation of the statistics where no correction for LD is attempted, in the style of [Pease & Hahn 2015].
#		The statistic within a window is simply (L - R)/(L + R) and the Z-score is simply (L - R)/sqrt(L + R).
#               The other option True used the blocks for block jackknifing that makes the estimator and its variance estimator robust to LD in the style of [Green 2010].
#		This block jackknifing  technique [Busin 1999] is described in a file called "wjack.pdf" by Nick Patterson you can find by Googling.
# output	The name of the file where results are printed, using one row per window if the parameter correction is False.
#		The default value is "" - then the results are not printed.
# progress	Whether the fuction prints progress to console or not. Defaults at True because I get restless otherwise.
#
# The function returns a list of lists of two numbers.
# The first one is the estimator of the Delta-statistic, the second is the Z-score i.e. the estimator divided by its estimated standard deviation.
def Delta(traw,
          populations,
          left,
          right,
          blocksize = 200000,
	  correction = False,
	  output = "",
	  progress = True):
    f = open(traw, "r")
    f.readline() # Removing the header.
    counts = [] # This will contain one list [L, R, N] per block.
    chromosome = -1 # Impossible value.
    L = 0
    R = 0
    N = 0
    counter = 0
    for line in f:
        if counter % 1000 == 0 and progress == True:
            print("Lines read: " + str(counter))
        counter = counter + 1
        if blocksize != 0: # If we're not using blocks at all, we just skip the next few lines.
            # If a block is ready, we save results and start a new block.
            if int(line.split()[0]) != chromosome:
                chromosome = int(line.split()[0])
                position = int(line.split()[3])
                counts.append([L, R, N])
                L = 0
                R = 0
                N = 0
            elif int(line.split()[3]) > position + blocksize:
                position = int(line.split()[3])
                counts.append([L, R, N])
                L = 0
                R = 0
                N = 0
        # Treat the variable.
        missing = False
        for word in left:
            product = 1
            for i in range(len(populations)):
                frequency = mymean([line.split()[j] for j in populations[i]])
                if frequency == "NA":
                    product = 0
                    missing = True
                else:
                    if word[i] == "B":
                        product = product*frequency
                    else:
                        product = product*(1 - frequency)
            L = L + product
        for word in right:
            product = 1
            for i in range(len(populations)):
                frequency = mymean([line.split()[j] for j in populations[i]])
                if frequency == "NA":
                    product = 0
                    missing = True
                else:
                    if word[i] == "B":
                        product = product*frequency
                    else:
                        product = product*(1 - frequency)
            R = R + product
        if missing == False:
            N = N + 1
    f.close()
    counts.append([L, R, N])
    if blocksize != 0:
        counts = counts[1:] # The first element is just zeroes so we remove it.
    # Compute the estimators.
    if correction == True: # Using the same block technique as ADMIXTOOLS. This technique should reduce biases caused by LD.
        zerodivision = False
        LminusR = sum([counts[i][0] for i in range(len(counts))]) - sum([counts[i][1] for i in range(len(counts))])
        LplusR = sum([counts[i][0] for i in range(len(counts))]) + sum([counts[i][1] for i in range(len(counts))])
        if LplusR > 0:
            Deltahat = LminusR/LplusR
        Deltahats = []
        for i in range(len(counts)):
            rang = [i for i in range(len(counts))]
            del rang[i]
            LminusR = sum([counts[i][0] for i in rang]) - sum([counts[i][1] for i in rang])
            LplusR = sum([counts[i][0] for i in rang]) + sum([counts[i][1] for i in rang])
            if LplusR > 0:
                Deltahats.append(LminusR/LplusR)
            else:
                zerodivision = True
        if zerodivision == False:
            # We should remove blocks that didn't have any observations.
            obs = []
            for i in range(len(counts)):
                if counts[i][2] > 0:
                    obs.append(i)
            counts = [counts[i] for i in obs]
            Deltahats = [Deltahats[i] for i in obs]
            # Then just plug into formulae.
            N = sum([counts[i][2] for i in range(len(counts))])
            DeltahatJ = 0
            taus = []
            for i in range(len(counts)):
                DeltahatJ = DeltahatJ + Deltahat + (counts[i][2]/N - 1)*Deltahats[i]
                taus.append(Deltahat*N/counts[i][2] - Deltahats[i]*(N/counts[i][2] - 1))
            # One more round!
            varhat = 0
            for i in range(len(counts)):
                varhat = varhat + ((taus[i] - DeltahatJ)**2)/(N/counts[i][2] - 1)
            varhat = varhat/len(counts)
            Zhat = DeltahatJ/math.sqrt(varhat)
        else:
            DeltahatJ = "NA"
            Zhat = "NA"
        result = [[DeltahatJ, Zhat]]
    else: # Not doing any corrections, just ignoring LD.
        result = []
        for i in range(len(counts)):
            if counts[i][0] + counts[i][1] > 0:
                result.append([(counts[i][0] - counts[i][1])/(counts[i][0] + counts[i][1]), (counts[i][0] - counts[i][1])/math.sqrt(counts[i][0] + counts[i][1])])
            else:
                result.append(["NA", "NA"])
    if output != "":
        g = open(output, "w")
        for window in result:
            nothing = g.write(str(window[0]) + " " + str(window[1]) + "\n")
    return(result)

########## SHORTCUTS ################################################################################################

# D
# The classic Patterson's D-statistic [Green 2010].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def D(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["ABBA", "BAAB"]
    right = ["BABA", "ABAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DFO
# The D_FO-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DFO(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "BBBAA", "ABABA", "AAABA", "ABABB", "AAABB", "BABAB", "BBBAB"]
    right = ["BAABA", "BBABA", "ABBAA", "AABAA", "ABBAB", "AABAB", "BAABB", "BBABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DIL
# The D_IL-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DIL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["ABBAA", "BBBAA", "BAABA", "AAABA", "BAABB", "AAABB", "ABBAB", "BBBAB"]
    right = ["ABABA", "BBABA", "BABAA", "AABAA", "BABAB", "AABAB", "ABABB", "BBABB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DFI
# The D_FI-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DFI(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "BABBA", "ABABA", "ABAAA", "ABABB", "ABAAB", "BABAB", "BABBB"]
    right = ["ABBAA", "ABBBA", "BAABA", "BAAAA", "BAABB", "BAAAB", "ABBAB", "ABBBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DOL
# The D_OL-statistic from the D_FOIL -method [Pease 2015].
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BAABA", "BABBA", "ABBAA", "ABAAA", "ABBAB", "ABAAB", "BAABB", "BABBB"]
    right = ["ABABA", "ABBBA", "BABAA", "BAAAA", "BABAB", "BAAAB", "ABABB", "ABBBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS1S8
# The ''freshman sum'' of Delta_S1 and minus Delta_S8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "AAABB", "BBBAA"]
    right = ["BAABA", "ABBAB", "AABAB", "BBABA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS2S8
# The ''freshman sum'' of Delta_S2 and minus Delta_S8.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["ABBAA", "BAABB", "AAABB", "BBBAA"]
    right = ["ABABA", "BABAB", "AABAB", "BBABA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS3S7
# The ''freshman sum'' of Delta_S3 and minus Delta_S7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABAAB", "BABBA"]
    right = ["ABBAA", "BAABB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS4S7
# The ''freshman sum'' of Delta_S4 and minus Delta_S7.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABAAB", "BABBA"]
    right = ["ABABA", "BABAB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS7S9
# The ''freshman sum'' of Delta_S7 and Delta_S9.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BAAAB", "ABBBA", "BAAAA", "ABBBB"]
    right = ["ABAAB", "BABBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS8S0
# The ''freshman sum'' of Delta_S8 and Delta_S10.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["AABAB", "BBABA", "AABAA", "BBABB"]
    right = ["AAABB", "BBBAA", "AAABA", "BBBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS3S4S7S9
# The ''freshman sum'' of Delta_S3, Delta_S4, minus Delta_S7 and Delta_S9.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "BAABA", "ABBAB", "ABAAB", "BABBA", "BAAAA", "ABBBB"]
    right = ["ABBAA", "BAABB", "ABABA", "BABAB", "BAAAB", "ABBBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DS1S2S8S0
# The ''freshman sum'' of Delta_S1, Delta_S2, minus Delta_S8 and Delta_S0.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABBAA", "BAABB", "AAABB", "BBBAA", "AABAA", "BBABB"]
    right = ["BAABA", "ABBAB", "ABABA", "BABAB", "AABAB", "BBABA", "AAABA", "BBBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA1A2
# The ''freshman sum'' of Delta_A1 and minus Delta_A2.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABABA", "BABAB"]
    right = ["ABBAA", "BAABB", "BAABA", "ABBAB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA1A3
# The ''freshman sum'' of Delta_A1 and minus Delta_A3.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "ABAAB", "BABBA"]
    right = ["ABBAA", "BAABB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA2A3
# The ''freshman sum'' of Delta_A2 and minus Delta_A3.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BAABA", "ABBAB", "ABAAB", "BABBA"]
    right = ["ABABA", "BABAB", "BAAAB", "ABBBA"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))

# DA1A2A3A4
# The ''freshman sum'' of Delta_A1, Delta_A2, minus Delta_A3 and Delta_A4.
# I'm not assuming the fourth population carries allele "A", that assumption can be made by filtering the data if desired.
def DOL(traw, populations, blocksize = 200000, correction = False, output = "", progress = True):
    left = ["BABAA", "ABABB", "BAABA", "ABBAB", "ABAAB", "BABBA", "BAAAA", "ABBBB"]
    right = ["ABBAA", "BAABB", "ABABA", "BABAB", "BAAAB", "ABBBA", "ABAAA", "BABBB"]
    return(Delta(traw, populations, left, right, blocksize, correction, output, progress))
