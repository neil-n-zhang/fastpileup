# ALLELECOUNT SCRIPT
# Written by Nathaniel Lim (May 2014) for Joshua Mell and Rosie Redfield.
# Written for Python 3
# Modified by Neil Zhang (Dec 2023)
# Output columns are: chr pos A C G T N a c g t n IN DEL indel
# Last field is only written for positions that indels were detected (i.e. leaves last column ragged).
# Indel field specifies counts of each insertion detected following the position and each deletion including the position.


# Module Imports
import sys

# ----------------------------------------------------------------------------------------------------
# Base Counting Subroutine *[Completed]
def base_counter(InputRow):
    InputList = []
    InputList = InputRow.split(sep = '\t')

    if InputList[2] in ["A","T","C","G"]:
        InputList[4]=InputList[4].replace(".",InputList[2]).replace(",",InputList[2].lower())
    
    # Cleaning up Base String + Indel Counting
    CleanString = ''
    countIn = 0
    countDel = 0
    IndelHolder = []
    IndelDeterminant = 0
    CleanBool = False
    
    for currentIndex, Strholder in enumerate(InputList[4]):
        # Skipping of '^' Signage
        if CleanBool == True:
            CleanBool = False
            continue
        
        if Strholder == '^':
            CleanBool = True
            continue
        
        # Skipping Indel
        if IndelDeterminant > 0:
            IndelDeterminant -= 1
            continue
        
        if Strholder == '+':
            countIn += 1
            
            # Determining Length of Indel
            # Since Illumina NGS has an upper limit of less than 999bp read length
            IndelDeterminant = 0
            
            if (currentIndex + 4) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 4].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 4]) + 3
                elif InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            elif (currentIndex + 3) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            
            IndelHolder.append(InputList[4][currentIndex:currentIndex + IndelDeterminant + 1])
            continue
        
        if Strholder == '-':
            countDel += 1
            
            # Determining Length of Indel
            # Since Illumina NGS has an upper limit of less than 999bp read length
            IndelDeterminant = 0
            
            if (currentIndex + 4) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 4].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 4]) + 3
                elif InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            elif (currentIndex + 3) <= len(InputList[4]):
                if InputList[4][currentIndex + 1: currentIndex + 3]. isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 3]) + 2
                elif InputList[4][currentIndex + 1: currentIndex + 2].isnumeric() == True:
                    IndelDeterminant = int(InputList[4][currentIndex + 1: currentIndex + 2]) + 1
            
            IndelHolder.append(InputList[4][currentIndex: currentIndex + len(str(IndelDeterminant)) + 1])
            continue
        
        CleanString += Strholder
    
    else:
    # Transferring Back Cleaned String
        InputList[4] = CleanString
        
        # '$' Signage Stripping
        InputList[4] = InputList[4].replace('$', '')
        
    # Base Count Var Initialization
    bigA = 0
    bigC = 0
    bigG = 0
    bigT = 0
    smallA = 0
    smallC = 0
    smallG = 0
    smallT = 0
    delBase = 0
    
    # Base Counting
    bigA = InputList[4].count('A')
    bigC = InputList[4].count('C')
    bigG = InputList[4].count('G')
    bigT = InputList[4].count('T')
    bigN = InputList[4].count('N')
    
    smallA = InputList[4].count('a')
    smallC = InputList[4].count('c')
    smallG = InputList[4].count('g')
    smallT = InputList[4].count('t')
    smallN = InputList[4].count('n')
    
    delBase = InputList[4].count('*')
    
    # Internal Check - Throws Out Error (NOT STD-IN/OUT Compatible: Should break pipeline)
    InternalCounter = 0
    InternalCounter = bigA + bigC + bigG + bigT + bigN + smallA + smallC + smallG + smallT + smallN + delBase
    if InternalCounter != int(InputList[3]):
        print('Error at position: ' + InputList[1])
        print('Reported count: ' + InputList[3])
        print('Internal counter: ' + str(InternalCounter))
        print('Internal sum: ' + str(len(InputList[4])))
        print('Number of insertions: ' + str(countIn))
        print('Number of deleions: ' + str(countDel))
        print('Post-processed bases: ' + InputList[4])
        sys.exit()
    
    # Indel Compilation
    IndelSetDict = set(IndelHolder)
    tmpIndelString = ''
    FinalIndelHolder = []
    
    for EveryIndel in IndelSetDict:
        tmpIndelString = ''
        tmpIndelString = str(IndelHolder.count(EveryIndel)) + ":" + EveryIndel
        FinalIndelHolder.append(tmpIndelString)
    
    # Return Output
    #FinalOutput = InputList[0] + '\t' + InputList[1] + '\t' + str(bigA) + '    ' + str(bigC) + '    ' + str(bigG) + '    ' + str(bigT) + '    ' + str(bigN) + '    ' + str(smallA) + '    ' + str(smallC) + '    ' + str(smallG) + '    ' + str(smallT) + '    ' + str(smallN) + '    ' + str(countIn) + '    ' + str(countDel) + '    ' + ';'.join(FinalIndelHolder)
    FinalOutput = "\t".join([InputList[0],InputList[1],str(bigA),str(bigC),str(bigG),str(bigT),str(bigN),str(smallA),str(smallC),str(smallG),str(smallT),str(smallN),str(countIn),str(countDel)])+ '\t' + ';'.join(FinalIndelHolder)
    
    return FinalOutput

def pileup2count(input_file,output_file):
    file_out=open(output_file,"w")
    with open(input_file) as file_in:
        for line in file_in:
            file_out.write(base_counter(line)+"\n")
    file_out.close()