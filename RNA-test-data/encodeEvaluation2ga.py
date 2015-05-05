#!/usr/bin/env python
from __future__ import division

import os
import sys
import errno
import string
import tarfile
import urllib2
import tempfile


def log(message):
    print >>sys.stderr, "== {}".format(message)


def get_tarball(tarball_url):
    tmp_tarball = tempfile.NamedTemporaryFile(mode='w+b')
    request = urllib2.urlopen(tarball_url)
    CHUNK = 16 * 1024
    log('Downloading tarball from {0}'.format(tarball_url))
    for chunk in iter(lambda: request.read(CHUNK), ''):
        if not chunk: break
        tmp_tarball.write(chunk)
        tmp_tarball.flush()
    return tmp_tarball


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
        # TODO: STAR does not report unique and multi split-maps, only total
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


def write_expression(analysisId, annotationId, quantfile, quantOutfile, tool='RSEM'):
    # RSEM gene expression table header:
    #   gene_id transcript_id(s)    length  effective_length    expected_count  TPM FPKM    pme_expected_count  
    #   pme_TPM pme_FPKM    TPM_ci_lower_bound  TPM_ci_upper_bound  FPKM_ci_lower_bound FPKM_ci_upper_bound
    # TODO: placeholder values need to be calculated then removed
    isNormalized = "True"
    units = "CPM"
    header = quantfile.readline()
    for expression in quantfile.readlines():
        fields = expression.split("\t")
        expressionLevel = fields[5]
        featureGroupId = fields[0]
        expressionId = analysisId
        rawCount = fields[4]
        score = (float(fields[10]) + float(fields[11]))/2
        outline = string.join([expressionId, annotationId, expressionLevel, featureGroupId, isNormalized, rawCount, str(score), units], "\t")
        quantOutfile.write("%s\n" % outline)


def write_rnaseq_tables(tar, description, annotationId):
    log("Writing rnaseq tables")
    for analysisId in get_analysis_ids_from_tar(tar):
        # create analysis id folder
        make_dir(analysisId)

        # output table
        rnaSeqTable = os.path.join(analysisId, "rnaseq.table")

        # write rnaseq table
        with open(rnaSeqTable, "w") as rnaQuantFile:
            writeRNAQuant(rnaQuantFile, analysisId, description, annotationId)


def write_counts_tables(tar):
    log("Writing counts tables")
    for analysisId, samstatsfile in get_members_from_tar(tar, 'Log.final.out'):
        # output table  
        countsTable = os.path.join(analysisId, "counts.table")

        samstats = getSamstats(samstatsfile)
        
        # write mapping stats table
        with open(countsTable, "w") as samOutfile:
            writeSamstats(samOutfile, samstats, analysisId)


def write_dist_tables(tar):
    log("Writing distribution tables")
    for analysisId, distfile in get_members_from_tar(tar, 'genome_cov'):
        # output table 
        distTable = os.path.join(analysisId, "dist.table")
        
        # write mapping distribution table
        distribution = getDistribution(distfile)
        with open(distTable, "w") as distOutfile:
            writeDistribution(distOutfile, distribution, analysisId, samstats["mapped"])


def write_expression_tables(tar, annotationId):
    log("Writing gene expression tables")
    for analysisId, quantfile in get_members_from_tar(tar, 'RSEM.gene'):
        # output table
        expTable = os.path.join(analysisId, "expression.table")
        
        # write expression table
        with open(expTable, "w") as quantOutfile:
            write_expression(analysisId, annotationId, quantfile, quantOutfile)


def main(argv):

    description = "RNAseq data from ENCODE evaluation"
    annotationId = "Gencodev16"

    tarball_url = 'http://genome.crg.es/~epalumbo/ENCODE-benchmark-data.tgz'
    tarball = get_tarball(tarball_url)
    
    tar = tarfile.open(tarball.name, mode='r')

    write_rnaseq_tables(tar, description, annotationId)
    write_counts_tables(tar)
    #TODO: get read distribution data
    #write_dist_tables(tar)
    write_expression_tables(tar, annotationId)

    tarball.close()
    tar.close()

    log("DONE")


if __name__ == '__main__':
    main(sys.argv)
