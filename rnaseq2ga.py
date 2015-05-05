from __future__ import division

import os
import sys
import errno
import string
import tarfile

def main(argv):

    # if len(argv) < 3:
    #     print 'usage: python %s description annotation' % sys.argv[0]
    #     sys.exit(1)

    description = "RNAseq data from ENCODE evaluation"
    annotationId = "Gencodev16"

    tar = tarfile.open('ENCODE-benchmark-data.tgz', 'r')

    #write_rnaseq_tables(tar, description, annotationId)
    write_counts_tables(tar)
    #write_dist_tables(tar)
    #write_expression_tables(tar, annotationId)

    tar.close()

def write_rnaseq_tables(tar, description, annotationId):
    for analysisId in get_analysis_ids_from_tar(tar):
        # create analysis id folder
        make_dir(analysisId)

        # output table
        rnaSeqTable = os.path.join(analysisId, "rnaseq.table")

        # write rnaseq table
        with open(rnaSeqTable, "w") as rnaQuantFile:
            writeRNAQuant(rnaQuantFile, analysisId, description, annotationId)


def write_counts_tables(tar):
    for analysisId, samstatsfile in get_members_from_tar(tar, 'Log.final.out'):
        # output table  
        countsTable = os.path.join(analysisId, "counts.table")

        samstats = getSamstats(samstatsfile)
        
        # write mapping stats table
        with open(countsTable, "w") as samOutfile:
            writeSamstats(samOutfile, samstats, analysisId)


def write_dist_tables(tar):
    for analysisId, distfile in get_members_from_tar(tar, 'genome_cov'):
        # output table 
        distTable = os.path.join(analysisId, "dist.table")
        
        # write mapping distribution table
        distribution = getDistribution(distfile)
        with open(distTable, "w") as distOutfile:
            writeDistribution(distOutfile, distribution, analysisId, samstats["mapped"])


def write_expression_tables(tar, annotationId):
    for analysisId, quantfile in get_members_from_tar(tar, 'RSEM.gene'):
        # output table
        expTable = os.path.join(analysisId, "expression.table")
        
        # write expression table
        with open(expTable, "w") as quantOutfile:
            write_expression(analysisId, annotationId, quantfile, quantOutfile)


def make_dir(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raise


def get_analysis_ids_from_tar(tar):
    for member in tar:
        if member.isdir():
            dirname, base = os.path.split(member.name)
            if not dirname:
                continue
            if os.path.split(dirname)[0]:
                continue
            yield base


def get_members_from_tar(tar, name):
    if not name:
        return
    for member in tar:
        if name in member.name:
            dirname = os.path.dirname(member.name)
            sample = os.path.basename(dirname)
            base = os.path.basename(member.name)
            member.name = "{0}_{1}".format(sample, base)
            yield (sample, tar.extractfile(member))


def write_expression(analysisId, annotationId, quantfile, quantOutfile, tool='RSEM'):
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

    def writecpms(analysisId, annotationId, isNormalized, units, cpmfile, cpmOutfile):
        #gene_id transcript_id(s)    length  effective_length    expected_count  TPM FPKM    pme_expected_count  
        #pme_TPM pme_FPKM    TPM_ci_lower_bound  TPM_ci_upper_bound  FPKM_ci_lower_bound FPKM_ci_upper_bound
        header = cpmfile.readline()
        for expression in cpmfile.readlines():
            fields = expression.split("\t")
            expressionLevel = fields[5]
            featureGroupId = fields[0]
            expressionId = analysisId
            rawCount = fields[4]
            score = (float(fields[10]) + float(fields[11]))/2
            outline = string.join([expressionId, annotationId, expressionLevel, featureGroupId, isNormalized, rawCount, str(score), units], "\t")
            cpmOutfile.write("%s\n" % outline)

    if tool == 'RSEM':
        #TODO: placeholder values need to be calculated then removed
        isNormalized = "True"
        units = "CPM"
        writecpms(analysisId, annotationId, isNormalized, units, quantfile, quantOutfile)
    elif tool == 'FLUX':
        #TODO: placeholder values need to be calculated then removed
        isNormalized = "True"
        units = "RPKM"
        writerpkms(analysisId, annotationId, isNormalized, units, quantfile, quantOutfile)
    elif tool == 'CUFFLINKS':
        #TODO: placeholder values need to be calculated then removed
        isNormalized = "True"
        units = "FPKM"
        writefpkms(analysisId, annotationId, isNormalized, units, quantfile, quantOutfile)


#TODO: placeholder values need to be calculated then removed
def getCount(expressionId):
    rawCount = 0
    return "%d" % rawCount

#TODO: placeholder values need to be calculated then removed
def getScore(expressionId):
    rawScore = 0.0

    return "%0.2f" % rawScore


def writeRNAQuant(outfile, analysisId, description, annotationId):
    """ Using placeholder value for the readGroupId
        Currently name is set to be same as analysisId
    """
    outline = string.join([analysisId, annotationId, description, analysisId, "readGroupId"], "\t")
    outfile.write("%s\n" % outline)


def getSamstats(file):
    content = {}
    for line in file.readlines():
        if '|' not in line:
            continue
        key, value = [item.strip() for item in line.split('|')]
        if key == "Number of input reads":
            content["readcount"] = value 
        elif key == 'Uniquely mapped reads number':    
            content["unique"] = value
        elif key in ['Number of reads mapped to multiple loci', 'Number of reads mapped to too many loci']:
            content["multi"] = content.get('multi', 0) + int(value)
        # STAR does not report unique and multi split-maps, only total
        elif key == 'Number of splices: Total':   
            content["msplice"] = content["usplice"] = value

    return content


def writeSamstats(outfile, contents, analysisId):
    outline = string.join([analysisId, str(contents["multi"]), contents["msplice"], contents["readcount"], contents["unique"], contents["usplice"]], "\t")
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
