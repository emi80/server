"""
Module responsible for handling protocol requests and returning
responses.
"""
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os
import json
import random

import ga4gh.protocol as protocol
import ga4gh.datamodel.references as references
import ga4gh.datamodel.reads as reads
import ga4gh.exceptions as exceptions
import ga4gh.datamodel.variants as variants
import ga4gh.datamodel.rna_quantification as rna_quantification


class AbstractBackend(object):
    """
    An abstract GA4GH backend.
    This class provides methods for all of the GA4GH protocol end points.
    """
    def __init__(self):
        self._variantSetIdMap = {}
        self._variantSetIds = []
        self._referenceSetIdMap = {}
        self._referenceSetIds = []
        self._readGroupSetIdMap = {}
        self._readGroupSetIds = []
        self._readGroupIds = []
        self._readGroupIdMap = {}
        self._rnaQuantificationIdMap = {}
        self._rnaQuantificationIds = []
        self._requestValidation = False
        self._responseValidation = False
        self._defaultPageSize = 100
        self._maxResponseLength = 2**20  # 1 MiB

    def getVariantSets(self):
        """
        Returns the list of VariantSets in this backend.
        """
        return list(self._variantSetIdMap.values())

    def getReadGroupSets(self):
        """
        Returns the list of ReadGroupSets in this backend.
        """
        return list(self._readGroupSetIdMap.values())

    def parsePageToken(self, pageToken, numValues):
        """
        Parses the specified pageToken and returns a list of the specified
        number of values. Page tokens are assumed to consist of a fixed
        number of integers seperated by colons. If the page token does
        not conform to this specification, raise a InvalidPageToken
        exception.
        """
        tokens = pageToken.split(":")
        # TODO define exceptions.InvalidPageToken and raise here.
        if len(tokens) != numValues:
            raise Exception("Invalid number of values in page token")
        # TODO catch a ValueError here when bad integers are passed and
        # convert this into the appropriate InvalidPageToken exception.
        values = map(int, tokens)
        return values

    def runSearchRequest(
            self, requestStr, requestClass, responseClass, objectGenerator):
        """
        Runs the specified request. The request is a string containing
        a JSON representation of an instance of the specified requestClass.
        We return a string representation of an instance of the specified
        responseClass in JSON format. Objects are filled into the page list
        using the specified object generator, which must return
        (object, nextPageToken) pairs, and be able to resume iteration from
        any point using the nextPageToken attribute of the request object.
        """
        self.startProfile()
        try:
            requestDict = json.loads(requestStr)
        except ValueError:
            raise exceptions.InvalidJsonException(requestStr)
        self.validateRequest(requestDict, requestClass)
        request = requestClass.fromJsonDict(requestDict)
        if request.pageSize is None:
            request.pageSize = self._defaultPageSize
        if request.pageSize <= 0:
            raise exceptions.BadPageSizeException(request.pageSize)
        responseBuilder = protocol.SearchResponseBuilder(
            responseClass, request.pageSize, self._maxResponseLength)
        nextPageToken = None
        for obj, nextPageToken in objectGenerator(request):
            responseBuilder.addValue(obj)
            if responseBuilder.isFull():
                break
        responseBuilder.setNextPageToken(nextPageToken)
        responseString = responseBuilder.getJsonString()
        self.validateResponse(responseString, responseClass)
        self.endProfile()
        return responseString

    def searchReadGroupSets(self, request):
        """
        Returns a GASearchReadGroupSetsResponse for the specified
        GASearchReadGroupSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadGroupSetsRequest,
            protocol.SearchReadGroupSetsResponse,
            self.readGroupSetsGenerator)

    def searchReads(self, request):
        """
        Returns a GASearchReadsResponse for the specified
        GASearchReadsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReadsRequest,
            protocol.SearchReadsResponse,
            self.readsGenerator)

    def searchReferenceSets(self, request):
        """
        Returns a GASearchReferenceSetsResponse for the specified
        GASearchReferenceSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferenceSetsRequest,
            protocol.SearchReferenceSetsResponse,
            self.referenceSetsGenerator)

    def searchReferences(self, request):
        """
        Returns a GASearchReferencesResponse for the specified
        GASearchReferencesRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchReferencesRequest,
            protocol.SearchReferencesResponse,
            self.referencesGenerator)

    def searchVariantSets(self, request):
        """
        Returns a GASearchVariantSetsResponse for the specified
        GASearchVariantSetsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantSetsRequest,
            protocol.SearchVariantSetsResponse,
            self.variantSetsGenerator)

    def searchVariants(self, request):
        """
        Returns a GASearchVariantsResponse for the specified
        GASearchVariantsRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchVariantsRequest,
            protocol.SearchVariantsResponse,
            self.variantsGenerator)

    def searchCallSets(self, request):
        """
        Returns a GASearchCallSetsResponse for the specified
        GASearchCallSetsRequest Object.
        """
        return self.runSearchRequest(
            request, protocol.SearchCallSetsRequest,
            protocol.SearchCallSetsResponse,
            self.callSetsGenerator)

    def searchRnaQuantification(self, request):
        """
        Returns a SearchRnaQuantificationResponse for the specified
        SearchRnaQuantificationRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchRnaQuantificationRequest,
            protocol.SearchRnaQuantificationResponse,
            self.rnaQuantificationGenerator)

    def searchExpressionLevel(self, request):
        """
        Returns a SearchExpressionLevelResponse for the specified
        SearchExpressionLevelRequest object.
        """
        return self.runSearchRequest(
            request, protocol.SearchExpressionLevelRequest,
            protocol.SearchExpressionLevelResponse,
            self.expressionLevelGenerator)

    # Iterators over the data hieararchy

    def _topLevelObjectGenerator(self, request, idMap, idList):
        """
        Generalisation of the code to iterate over the objects at the top
        of the data hierarchy.
        """
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = self.parsePageToken(request.pageToken, 1)
        while currentIndex < len(idList):
            objectId = idList[currentIndex]
            object_ = idMap[objectId]
            currentIndex += 1
            nextPageToken = None
            if currentIndex < len(idList):
                nextPageToken = str(currentIndex)
            yield object_.toProtocolElement(), nextPageToken

    def readGroupSetsGenerator(self, request):
        """
        Returns a generator over the (readGroupSet, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._readGroupSetIdMap, self._readGroupSetIds)

    def referenceSetsGenerator(self, request):
        """
        Returns a generator over the (referenceSet, nextPageToken) pairs
        defined by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._referenceSetIdMap, self._referenceSetIds)

    def variantSetsGenerator(self, request):
        """
        Returns a generator over the (variantSet, nextPageToken) pairs defined
        by the specified request.
        """
        return self._topLevelObjectGenerator(
            request, self._variantSetIdMap, self._variantSetIds)

    def readsGenerator(self, request):
        # Local utility functions to save some typing
        def getPosition(readAlignment):
            return readAlignment.alignment.position.position

        def getEndPosition(readAlignment):
            return getPosition(readAlignment) + \
                len(readAlignment.alignedSequence)

        if len(request.readGroupIds) != 1:
            raise exceptions.NotImplementedException(
                "Read search over multiple readGroups not supported")
        readGroupId = request.readGroupIds[0]
        try:
            readGroup = self._readGroupIdMap[request.readGroupIds[0]]
        except KeyError:
            raise exceptions.ReadGroupNotFoundException(readGroupId)
        startPosition = request.start
        equalPositionsToSkip = 0
        if request.pageToken is not None:
            startPosition, equalPositionsToSkip = self.parsePageToken(
                request.pageToken, 2)
        iterator = readGroup.getReadAlignments(
            request.referenceId, startPosition,
            request.end)
        readAlignment = next(iterator, None)
        if request.pageToken is not None:
            # First, skip any reads with position < startPosition
            # or endPosition < request.start
            while (getPosition(readAlignment) < startPosition or
                   getEndPosition(readAlignment) < request.start):
                readAlignment = next(iterator, None)
                if readAlignment is None:
                    raise exceptions.BadPageTokenException(
                        "Inconsistent page token provided")
            # Now, skip equalPositionsToSkip records which have position
            # == startPosition
            equalPositionsSkipped = 0
            while equalPositionsSkipped < equalPositionsToSkip:
                if getPosition(readAlignment) != startPosition:
                    raise exceptions.BadPageTokenException(
                        "Inconsistent page token provided")
                equalPositionsSkipped += 1
                readAlignment = next(iterator, None)
                if readAlignment is None:
                    raise exceptions.BadPageTokenException(
                        "Inconsistent page token provided")
        # The iterator is now positioned at the correct place.
        while readAlignment is not None:
            nextReadAlignment = next(iterator, None)
            nextPageToken = None
            if nextReadAlignment is not None:
                if getPosition(readAlignment) == \
                        getPosition(nextReadAlignment):
                    equalPositionsToSkip += 1
                else:
                    equalPositionsToSkip = 0
                nextPageToken = "{}:{}".format(
                    getPosition(nextReadAlignment), equalPositionsToSkip)
            yield readAlignment, nextPageToken
            readAlignment = nextReadAlignment

    def variantsGenerator(self, request):
        """
        Returns a generator over the (variant, nextPageToken) pairs defined by
        the specified request.
        """
        # TODO this method should also use the interval search semantics
        # in the readsGenerator above.
        if len(request.variantSetIds) != 1:
            raise exceptions.NotImplementedException(
                "VariantSearch search over multiple variantSets not supported")
        variantSetId = request.variantSetIds[0]
        try:
            variantSet = self._variantSetIdMap[request.variantSetIds[0]]
        except KeyError:
            raise exceptions.VariantSetNotFoundException(variantSetId)
        startPosition = request.start
        if request.pageToken is not None:
            startPosition, = self.parsePageToken(request.pageToken, 1)
        iterator = variantSet.getVariants(
            request.referenceName, startPosition, request.end,
            request.variantName, request.callSetIds)
        variant = next(iterator, None)
        while variant is not None:
            nextVariant = next(iterator, None)
            nextPageToken = None
            if nextVariant is not None:
                nextPageToken = "{}".format(nextVariant.start)
            yield variant, nextPageToken
            variant = nextVariant

    def callSetsGenerator(self, request):
        """
        Returns a generator over the (callSet, nextPageToken) pairs defined by
        the specified request.
        """
        if len(request.variantSetIds) != 1:
            raise exceptions.NotImplementedException(
                "Searching over multiple variantSets is not supported.")
        if request.name is not None:
            raise exceptions.NotImplementedException(
                "Searching over names is not supported")
        variantSetId = request.variantSetIds[0]
        try:
            variantSet = self._variantSetIdMap[variantSetId]
        except KeyError:
            raise exceptions.VariantSetNotFound(variantSetId)
        return self._topLevelObjectGenerator(
            request, variantSet.getCallSetIdMap(),
            variantSet.getCallSetIds())

    def rnaQuantificationGenerator(self, request):
        """
        Returns a generator over the (rnaQuantification, nextPageToken) pairs defined
        by the specified request.
        """
        rnaQuantificationId = request.rnaQuantificationId
        try:
            rnaQuant = self._rnaQuantificationIdMap[rnaQuantificationId]
        except KeyError:
            raise exceptions.RnaQuantificationNotFoundException(rnaQuantificationId)
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = self.parsePageToken(request.pageToken, 1)
        rnaQuantIterator = rnaQuant.getRnaQuantification(rnaQuantificationId)
        rnaQuantData = next(rnaQuantIterator, None)
        while rnaQuantData is not None:
            nextRnaQuantData = next(rnaQuantIterator, None)
            nextPageToken = None
            if nextRnaQuantData is not None:
                currentIndex += 1
                nextPageToken = "{}".format(currentIndex)
            rnaQuantification = protocol.RnaQuantification()
            rnaQuantification.annotationIds = rnaQuantData.annotationIds
            rnaQuantification.description = rnaQuantData.description
            rnaQuantification.id = rnaQuantData.id
            rnaQuantification.name = rnaQuantData.name
            rnaQuantification.readGroupId = rnaQuantData.readGroupId
            yield rnaQuantification, nextPageToken
            rnaQuantData = nextRnaQuantData

    def expressionLevelGenerator(self, request):
        expressionLevelId = request.expressionLevelId
        featureGroupId = request.featureGroupId
        rnaQuantificationId = request.rnaQuantificationId
        currentIndex = 0
        if request.pageToken is not None:
            currentIndex, = self.parsePageToken(request.pageToken, 1)
        if rnaQuantificationId is not None:
            rnaQuantificationIds = [rnaQuantificationId]
        else:
            rnaQuantificationIds = self._rnaQuantificationIds
        for rnaQuantId in rnaQuantificationIds:
            rnaQuant = self._rnaQuantificationIdMap[rnaQuantId]
            expressionLevelIterator = rnaQuant.getExpressionLevel(expressionLevelId, featureGroupId)
            expressionLevelData = next(expressionLevelIterator, None)
            while expressionLevelData is not None:
                nextExpressionLevelData = next(expressionLevelIterator, None)
                nextPageToken = None
                if nextExpressionLevelData is not None:
                    currentIndex += 1
                    nextPageToken = "{}".format(currentIndex)
                expressionLevel = protocol.ExpressionLevel()
                expressionLevel.annotationId = expressionLevelData.annotationId
                expressionLevel.expression = expressionLevelData.expression
                expressionLevel.featureGroupId = expressionLevelData.featureGroupId
                expressionLevel.id = expressionLevelData.id
                expressionLevel.isNormalized = expressionLevelData.isNormalized
                expressionLevel.rawReadCount = expressionLevelData.rawReadCount
                expressionLevel.score = expressionLevelData.score
                expressionLevel.units = expressionLevelData.units
                yield expressionLevel, nextPageToken
                expressionLevelData = nextExpressionLevelData

    def startProfile(self):
        """
        Profiling hook. Called at the start of the runSearchRequest method
        and allows for detailed profiling of search performance.
        """
        pass

    def endProfile(self):
        """
        Profiling hook. Called at the end of the runSearchRequest method.
        """
        pass

    def validateRequest(self, jsonDict, requestClass):
        """
        Ensures the jsonDict corresponds to a valid instance of requestClass
        Throws an error if the data is invalid
        """
        if self._requestValidation:
            if not requestClass.validate(jsonDict):
                raise exceptions.RequestValidationFailureException(
                    jsonDict, requestClass)

    def validateResponse(self, jsonString, responseClass):
        """
        Ensures the jsonDict corresponds to a valid instance of responseClass
        Throws an error if the data is invalid
        """
        if self._responseValidation:
            jsonDict = json.loads(jsonString)
            if not responseClass.validate(jsonDict):
                raise exceptions.ResponseValidationFailureException(
                    jsonDict, responseClass)

    def setRequestValidation(self, requestValidation):
        """
        Set enabling request validation
        """
        self._requestValidation = requestValidation

    def setResponseValidation(self, responseValidation):
        """
        Set enabling response validation
        """
        self._responseValidation = responseValidation

    def setDefaultPageSize(self, defaultPageSize):
        """
        Sets the default page size for request to the specified value.
        """
        self._defaultPageSize = defaultPageSize

    def setMaxResponseLength(self, maxResponseLength):
        """
        Sets the approximate maximum response length to the specified
        value.
        """
        self._maxResponseLength = maxResponseLength


class EmptyBackend(AbstractBackend):
    """
    A GA4GH backend that contains no data.
    """


class SimulatedBackend(AbstractBackend):
    """
    A GA4GH backend backed by no data; used mostly for testing
    """
    def __init__(self, randomSeed=0, numCalls=1, variantDensity=0.5,
                 numVariantSets=1):
        super(SimulatedBackend, self).__init__()
        self._randomSeed = randomSeed
        self._randomGenerator = random.Random()
        self._randomGenerator.seed(self._randomSeed)
        for i in range(numVariantSets):
            variantSetId = "simVs{}".format(i)
            seed = self._randomGenerator.randint(0, 2**32 - 1)
            variantSet = variants.SimulatedVariantSet(
                seed, numCalls, variantDensity, variantSetId)
            self._variantSetIdMap[variantSetId] = variantSet
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # Reads
        readGroupSetId = "aReadGroupSet"
        readGroupSet = reads.SimulatedReadGroupSet(readGroupSetId)
        self._readGroupSetIdMap[readGroupSetId] = readGroupSet
        for readGroup in readGroupSet.getReadGroups():
            self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())


class FileSystemBackend(AbstractBackend):
    """
    A GA4GH backend backed by data on the file system
    """
    def __init__(self, dataDir):
        super(FileSystemBackend, self).__init__()
        self._dataDir = dataDir
        # TODO this code is very ugly and should be regarded as a temporary
        # stop-gap until we deal with iterating over the data tree properly.
        # Variants
        variantSetDir = os.path.join(self._dataDir, "variants")
        for variantSetId in os.listdir(variantSetDir):
            relativePath = os.path.join(variantSetDir, variantSetId)
            if os.path.isdir(relativePath):
                self._variantSetIdMap[variantSetId] = \
                    variants.HtslibVariantSet(variantSetId, relativePath)
        self._variantSetIds = sorted(self._variantSetIdMap.keys())

        # References
        referenceSetDir = os.path.join(self._dataDir, "references")
        for referenceSetId in os.listdir(referenceSetDir):
            relativePath = os.path.join(referenceSetDir, referenceSetId)
            if os.path.isdir(relativePath):
                referenceSet = references.ReferenceSet(
                    referenceSetId, relativePath)
                self._referenceSetIdMap[referenceSetId] = referenceSet
        self._referenceSetIds = sorted(self._referenceSetIdMap.keys())

        # Reads
        readGroupSetDir = os.path.join(self._dataDir, "reads")
        for readGroupSetId in os.listdir(readGroupSetDir):
            relativePath = os.path.join(readGroupSetDir, readGroupSetId)
            if os.path.isdir(relativePath):
                readGroupSet = reads.HtslibReadGroupSet(
                    readGroupSetId, relativePath)
                self._readGroupSetIdMap[readGroupSetId] = readGroupSet
                for readGroup in readGroupSet.getReadGroups():
                    self._readGroupIdMap[readGroup.getId()] = readGroup
        self._readGroupSetIds = sorted(self._readGroupSetIdMap.keys())
        self._readGroupIds = sorted(self._readGroupIdMap.keys())

        #Rna Quantification
        rnaQuantDir = os.path.join(self._dataDir, "rnaQuant")
        for rnaQuantId in os.listdir(rnaQuantDir):
            relativePath = os.path.join(rnaQuantDir, rnaQuantId)
            if os.path.isdir(relativePath):
                rnaQuantification = rna_quantification.RNASeqResult(
                    rnaQuantId, relativePath)
                self._rnaQuantificationIdMap[rnaQuantId] = rnaQuantification
        self._rnaQuantificationIds = sorted(self._rnaQuantificationIdMap.keys())

