
class TF1(object):
    PValue=0
    SNPName=1
    SNPChr=2
    SNPChrPos=3
    ProbeName=4
    ProbeChr=5
    ProbeCenterChrPos=6
    CisTrans=7
    SNPType=8
    AlleleAssessed=9
    OverallZScore=10
    DatasetsWhereSNPProbePairIsAvailableAndPassesQC=11
    DatasetsZScores=12
    HUGO=13
    FDR=14

def alleles(TF, comps):
    alleles = {a for a in comps[TF.SNPType].split("/")}

    #TODO: figure this out
    assessed = comps[TF.AlleleAssessed]
    for a in alleles:
        if a == assessed:
            effect_allele = a
        else:
            reference_allele = a
    return reference_allele, effect_allele
