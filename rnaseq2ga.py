import os
import sys
import string
import tarfile

def main(argv):

    if len(argv) < 3:
        print 'usage: python %s rnaQuantId description annotation' % sys.argv[0]
        sys.exit(1)

    analysisIds = argv[1].split(',')
    description = argv[2]
    annotationId = argv[3]

    for analysisId in analysisIds:
        if not analysisId:
            continue
        try:
            os.mkdir(analysisId)
        except:
            pass
        # input files
        mapStats = "{0}.mappedReads.txt".format(analysisId)
        genomeCov = "{0}.genome.covg.tsv".format(analysisId)
        quantTable = "{0}.transcript.gtf".format(analysisId)

        # output tables
        rnaSeqTable = os.path.join(analysisId, "rnaseq.table")
        countsTable = os.path.join(analysisId, "counts.table")
        distTable = os.path.join(analysisId, "dist.table")
        expTable = os.path.join(analysisId, "expression.table")

        # write rnaseq table
        with open(rnaSeqTable, "w") as rnaQuantFile:
            writeRNAQuant(rnaQuantFile, analysisId, description)

        # write mapping stats table
        with open(mapStats, "r") as samstatsfile:
            samstats = getSamstats(samstatsfile)
        with open(countsTable, "w") as samOutfile:
            writeSamstats(samOutfile, samstats, analysisId)

        # write mapping distribution table
        with open(genomeCov, "r") as distfile:
            distribution = getDistribution(distfile)
        with open(distTable, "w") as distOutfile:
            writeDistribution(distOutfile, distribution, analysisId, samstats["mapped"])

        # write expression table
        #TODO: placeholder values need to be calculated then removed
        isNormalized = "True"
        units = "RPKM"
        with open(quantTable, "r") as quantfile, open(expTable, "w") as quantOutfile:
            writerpkms(analysisId, annotationId, isNormalized, units, quantfile, quantOutfile)

def writerpkms(analysisId, annotationId, isNormalized, units, rpkmfile, rpkmOutfile):
    for transcript in rpkmfile.readlines():
        fields = transcript.split("\t")
        content = {}
        content['seqname'] = fields[0]
        content['source'] = fields[1]
        content['feature'] = fields[2]
        content['start'] = fields[3]
        content['end'] = fields[4]
        content['score'] = fields[5]
        content['strand'] = fields[6]
        content['frame'] = fields[7]
        attrList = (attr.strip() for attr in filter(None, fields[8].strip().split(';')))
        content['attribute'] = dict((attr.split()[0].strip(), attr.split()[1].strip("\"")) for attr in attrList)

        expressionLevel = content['attribute']['RPKM']
        featureGroupId =  content['attribute']['transcript_id']
        expressionId = analysisId + "-{seqname}:{start}-{end}".format(**content)
        rawCount = content['attribute']['reads']
        score = content['score']
        outline = string.join([expressionId, annotationId, expressionLevel, featureGroupId, isNormalized, rawCount, score, units], "\t")
        rpkmOutfile.write("%s\n" % outline)

def writefpkms(analysisId, annotationId, isNormalized, units, fpkmfile, fpkmOutfile):
    header = fpkmfile.readline()
    for expression in fpkmfile.readlines():
        fields = expression.split("\t")
        expressionLevel = fields[9]
        featureGroupId = fields[4]
        expressionId = fields[0]
        rawCount = getCount(expressionId)
        score = getScore(expressionId)
        outline = string.join([expressionId, annotationId, expressionLevel, featureGroupId, isNormalized, rawCount, score, units], "\t")
        fpkmOutfile.write("%s\n" % outline)

#TODO: placeholder values need to be calculated then removed
def getCount(expressionId):
    rawCount = 0
    return "%d" % rawCount

#TODO: placeholder values need to be calculated then removed
def getScore(expressionId):
    rawScore = 0.0

    return "%0.2f" % rawScore


def writeRNAQuant(outfile, analysisId, description):
    """ Using placeholder values for the annotations and readGroupId
        Currently name is set to be same as analysisId
    """
    annotations = string.join(["1", "2"], ",")
    outline = string.join([analysisId, annotations, description, analysisId, "readGroupId"], "\t")
    outfile.write("%s\n" % outline)


def getSamstats(file):
    content = {}
    _ = file.readline()
    content["total"] = file.readline().split(":")[1].strip()
    content["mapped"] = file.readline().split(":")[1].strip()
    content["pairs"] = file.readline().split(":")[1].strip()

    return content


def writeSamstats(outfile, contents, analysisId):
    #outline = string.join([analysisId, contents["multi"], contents["msplice"], contents["readcount"], contents["unique"], contents["usplice"]], "\t")
    outline = string.join([analysisId, contents["mapped"], contents["mapped"], contents["total"], contents["mapped"], contents["mapped"]], "\t")
    outfile.write("%s\n" % outline)


def getDistribution(file):
    content = {}
    for line in file.readlines():
        split = line.split('\t')
        if split[3] == 'total':
            continue
        content[split[2]] = content.get(split[2],0) + int(split[3])
    return content


def writeDistribution(outfile, contents, analysisId, fraction):
    outline = string.join([analysisId, str(contents["exon"]), fraction, str(contents["intergenic"]), str(contents["intron"])], "\t")
    outfile.write("%s\n" % outline)


if __name__ == '__main__':
    main(sys.argv)
