#!/usr/bin/env python

import networkx as nx
import re, math, os, sys, operator, random
from scipy.sparse import coo_matrix
# from scipy.linalg import expm
import imp
#import kernel_scipy
from optparse import OptionParser
import csv
import numpy as np
#schema=imp.load_source('ndexSchema','ndexSchema.py')
import ndexSchema as schema
#util=imp.load_source('ndexUtil','ndexUtil.py')
import ndexUtil as util
#nc=imp.load_source('ndexClient','ndexClient.py')
import ndex.client as nc
#kernel=imp.load_source('kernel_scipy','kernel_scipy.py')
import kernel_scipy as kernel
import operator

class TableEntityMapper:
    """This class is a container for a mapping table between two namespaces.
    Expects a filename for a file that is two tab delimited columns, source\ttarget
    Retutns a list (!) of ids"""

    def __init__(self, filename, missingStrList=[]):
        self.mappingTable = {}
        self.readFileToMappingTable(filename)
        self.missingStrList = missingStrList

    def readFileToMappingTable(self, filename):
        with open(filename, "r") as ins:
            for line in ins:
                arr = line.split('\t')
                if arr[0] not in self.mappingTable:
                    self.mappingTable[arr[0]] = [arr[1].rstrip()]
                else:
                    print("multiple mappings of %s detected" % arr[0])
                    self.mappingTable[arr[0]].append(arr[1])


    def map(self, entity):
        try:
            return self.mappingTable[entity]
        except:
            return self.missingStrList


def mapByPrefix(i, MGImapper=None, RGDmapper=None, HGNCmapper=None, targetPrefix='hgnc:'):
    st = i.lower()
    if st.startswith(targetPrefix):
        hgout = [i[len(targetPrefix):len(i)]]
    elif st.startswith('hgnc:') and HGNCmapper is not None:
        hgout = HGNCmapper.map(i[5:len(i)])
    elif st.startswith('mgi:') and MGImapper is not None:
        hgout = MGImapper.map(i[4:len(i)])
    elif st.startswith('rgd:') and RGDmapper is not None:
        hgout = RGDmapper.map(i[4:len(i)])
    else:
        hgout = [i]
    return hgout


def filterByPrefix(l, targetPrefix='hgnc:'):
    out = []
    for i in l:
        st = i.lower()
        if not st.startswith('bel:'):
            if st.startswith(targetPrefix):
                out.append(i[len(targetPrefix):len(i)])
            else:
                out.append(i)
    return out


class NdexToGeneSif(util.NetworkWrapper):
    def __init__(self, ndexNetwork, prefix='hgnc:', MGImapper=None, RGDmapper=None, HGNCmapper=None):
        self.prefix = prefix
        self.network = ndexNetwork
        self.supportToEdgeMap = {}
        self.citationToSupportMap = {}
        self.nodeLabelMap = {}
        self.termLabelMap = {}
        #new vars
        self.basetermMap = {}
        self.missingBaseterms = []
        #fix this guy for non-bel:
        self.nodeBaseterms = {}
        self.MGImapper = MGImapper
        self.RGDmapper = RGDmapper
        self.HGNCmapper = HGNCmapper

        for nodeId, node in ndexNetwork['nodes'].iteritems():
            self.nodeLabelMap[int(nodeId)] = self.getNodeLabel(node)
            self.nodeBaseterms[int(nodeId)] = self.getNodeBaseterms(node)

        for edge in ndexNetwork['edges'].values():
            for supportId in edge['supportIds']:
                supports = ndexNetwork['supports']
                support = supports[str(supportId)]
                if supportId in self.supportToEdgeMap:
                    edgeList = self.supportToEdgeMap[supportId]
                else:
                    edgeList = []
                edgeList.append(edge)
                self.supportToEdgeMap[supportId] = edgeList

        for supportId in self.supportToEdgeMap.keys():
            support = ndexNetwork['supports'][str(supportId)]
            citationId = support['citationId']
            if citationId in self.citationToSupportMap:
                supportIdList = self.citationToSupportMap[citationId]
            else:
                supportIdList = []
            supportIdList.append(supportId)
            self.citationToSupportMap[citationId] = supportIdList

    #this is new
    def edgeExpandToSIF(self, edge):
        """Takes an edge and simplifies to gene symbol and interaction type.
       All baseterms associated with subject inherit all interactions with all objects"""
        subjectBasetermList = []
        objectBasetermList = []
        predicateLabel = "missing"
        subjectId = edge['subjectId']
        objectId = edge['objectId']
        if subjectId in self.nodeBaseterms:
            subjectBasetermList=self.nodeBaseterms[subjectId]
        if objectId in self.nodeBaseterms:
            objectBasetermList=self.nodeBaseterms[objectId]

        predicateId = edge['predicateId']
        predicateLabel = util.stripPrefixes(self.getTermLabel(predicateId), targetPrefix=self.prefix)
        sifTriples = []

        subjMapped = []
        for subj in subjectBasetermList:
            subjMapped.extend(mapByPrefix(subj, self.MGImapper, self.RGDmapper, self.HGNCmapper, targetPrefix=self.prefix))
        subjMapped = filterByPrefix(subjMapped)

        objMapped = []
        for obj in objectBasetermList:
            objMapped.extend(mapByPrefix(obj, self.MGImapper, self.RGDmapper, self.HGNCmapper, targetPrefix=self.prefix))
        objMapped = filterByPrefix(objMapped)

        for subj in subjMapped:
            for obj in objMapped:
                sifTriples.append([subj, predicateLabel, obj])

        return sifTriples


    def getNodeBaseterms(self, node):
        if 'name' in node and node['name']:
            return [node['name']]

        elif 'represents' in node:
            return self.basetermMap[node['represents']]

        else:
            return ["node %s"] % (node['id'])


    def writeSIF(self, fileName=None, append=False):
        if fileName and not append:
            output = open(fileName, 'w')
        elif fileName and append:
            output = open(fileName, 'a')
        else:
            output = sys.stdout

        for edge in self.network['edges'].values():
            sifs = self.edgeExpandToSIF(edge)
            for s in sifs:
                if len(s) == 3:
                    output.write('\t'.join(s) + '\n')
                
        if fileName:
            output.close

    def edgeList(self):
        edge_list=list()

        for edge in self.network['edges'].values():
            sifs = self.edgeExpandToSIF(edge)
            for s in sifs:
                edge_list.append((s[0],s[2]))
        return edge_list

    def checkType(self,term):
        if 'name' in term.keys():
            return 'baseterm'
        elif 'functionTermId' in term.keys():
            return 'functionterm'
        elif 'edgeId' in term.keys():
            return 'reifiededgeterm'
        else:
            return 'unknown'

    def getTermLabel(self, termId):
        if termId in self.termLabelMap:
            return self.termLabelMap[termId]
        else:
            label = "error"
            term = self.getTermById(termId)
            #type = term['type'].lower()
            type=self.checkType(term)
            if type == "baseterm":
                name = term['name']
                if 'namespaceId' in term and term['namespaceId']:
                    namespaceId = term['namespaceId']
                    namespace = self.getNamespaceById(namespaceId)

                    if namespace:
                        if namespace['prefix']:
                            label = "%s:%s" % (namespace['prefix'], name)
                        elif namespace['uri']:
                            label = "%s%s" % (namespace['uri'], name)
                        else:
                            label = "%s%s" % ('BEL:', name)
                    else:
                        label = "%s%s" % ('BEL:', name)
                else:
                    label = "%s%s" % ('BEL:', name)
                self.basetermMap[termId] = [label]

            elif type == "functionterm":
                functionTermId = term['functionTermId']
                functionLabel = self.getTermLabel(functionTermId)
                functionLabel = util.getFunctionAbbreviation(functionLabel)
                parameterLabels = []
                if (termId not in self.basetermMap):
                    self.basetermMap[termId] = []
                for parameterId in term['parameterIds']:
                    parameterLabel = self.getTermLabel(parameterId)
                    parameterLabels.append(parameterLabel)
                    if parameterId in self.basetermMap:
                        self.basetermMap[termId].extend(self.basetermMap[parameterId])
                    else:
                        self.missingBaseterms.append(parameterId)
                label = "%s(%s)" % (functionLabel, ",".join(parameterLabels))

            elif type == "reifiededgeterm":
                edgeId = term['edgeId']
                edges = self.network['edges']
                self.basetermMap[termId] = []
                if edgeId in edges:
                    reifiedEdge = edges[edgeId]
                    label = "(%s)" % (self.getEdgeLabel(reifiedEdge))
                else:
                    label = "(reifiedEdge: %s)" % (edgeId)

            else:
                label = "term: %s" % (termId)
                self.basetermMap[termId] = []

            self.termLabelMap[termId] = label
            return label


def queryVector(query, labels):
    out = {}
    weight_sum = 0.0
    for i in labels:
        if i in query:
            out[i] = 1.0
            weight_sum += 1.0
        else:
            out[i] = 0.0

    for i in labels:
        out[i] = out[i] / weight_sum

    return out


def weightedQueryVector(query, labels, weights):
    out = {}
    weight_sum = 0
    for i in labels:
        if i in query:
            try:
                out[i] = weights[i]
                weight_sum += weights[i]
            except KeyError:
                out[i] = 1
                weight_sum += 1
        else:
            out[i] = 0

    for i in labels:
        out[i] = np.float32(out[i] / weight_sum)

    return out


def readNodeWeights(filename):
    weights = {}
    for l in open(filename, 'r'):
        line=l.rstrip().split('\t')
        if len(line)==2:
            weights[line[0]] = float(line[1])
        else:
            print "problem line:"+'TAB'.join(line)
    return weights


def readSif(filename):
    triples = []
    for line in csv.reader(open(filename, 'rb'), delimiter='\t'):
        triples.append(line)
    return triples


def filterSif(triples, scores, desired_nodes=30):
    scores_made_the_cut = sorted(scores.items(), key=operator.itemgetter(1), reverse=True)[0:(desired_nodes)]
    nodes_made_the_cut = [i[0] for i in scores_made_the_cut]
    edges_made_the_cut = []
    for tr in triples:
        if (tr[0] in nodes_made_the_cut) and (tr[2] in nodes_made_the_cut):
            edges_made_the_cut.append(tr)
    return edges_made_the_cut


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-q", "--query", dest="query", action="store", type="string", default=None,
                      help="File containining a list of query genes")
    parser.add_option("-s", "--sif", dest="sif", action="store", type="string", default="out.sif",
                      help="Output file for intermediate sif, translated to HGNC")
    parser.add_option("-k", "--kernel", dest="kernel", action="store", type="string", default="kernel.txt",
                      help="Output file for kernel")
    parser.add_option("-t", "--diffusion-time", dest="diffusion_time", type="float", default=0.1,
                      help="Time parameter (scale) for diffusion.  Default is 0.1")
    parser.add_option("-d", "--diffused-query", dest="diffused_query", type="string", default="diffused.txt",
                      help="Output file for diffused query.")
    parser.add_option("-u", "--username", dest="username", type="string", default=None, help="NDEx username")
    parser.add_option("-p", "--password", dest="password", type="string", default=None, help="NDEx password")
    parser.add_option("-n", "--network-uuid", dest="network_uuid", type="string",
                      default="9ea3c170-01ad-11e5-ac0f-000c29cb28fb",
                      help="uuid of query, default is for the BEL large corpus")
    parser.add_option("-m", "--mgi-mapping-table", dest="mgi", type="string", default="bel_MGItoHGNC.tab",
                      help="file with MGI to target species mapping table")
    parser.add_option("-r", "--rgd-mapping-table", dest="rgd", type="string", default="bel_RGDtoHGNC.tab",
                      help="file with RGC to target species mapping table")
    parser.add_option("-H", "--hgnc-mapping-table", dest="hgnc", type="string", default=None,
                      help="file with HGNC to target species mapping table")
    parser.add_option("-b", "--bel-file", dest="bel", type="string", default=None)
    parser.add_option("-S", "--species", dest="species", type="string", default="human",
                      help="Target species to map to.  Default is 'human', other options right now are 'rat' and 'mouse'")
    parser.add_option("-x", "--do-not-map", dest="nomap", action="store_true", default=False,
                      help="Do not do orthology mappings")
    parser.add_option("-w", "--weighted-query", dest="weighted_query", default=None,
                      help = "weighted query file.  Expects a column of IDs and a column of scores, no header")
    parser.add_option("-N", "--number-of-entities", dest="entities", type="int", default=30, help="Number of entities to have in the final sif, default  is 30")
    parser.add_option("-o", "--out-sif", dest="filtered_sif", type="string", default="filtered.sif", help="Filtered output sif.")
    (opts, args) = parser.parse_args()


    #lib_path = os.path.abspath('/Users/danielcarlin/projects/TieDIE/lib')
    #sys.path.append(lib_path)

    #establishes orthology mappers
    if (opts.nomap):
        MGImapper = None
        RGDmapper = None
        HGNCmapper = None
        prefix = "NO_LEGAL_PREFIX"
    elif (opts.species == 'human'):
        MGImapper = TableEntityMapper(opts.mgi)
        RGDmapper = TableEntityMapper(opts.rgd)
        HGNCmapper = None
        prefix = 'hgnc:'
    elif (opts.species == 'mouse'):
        HGNCmapper = TableEntityMapper(opts.hgnc)
        RGDmapper = TableEntityMapper(opts.rgd)
        MGImapper = None
        prefix = 'mgi:'
    elif (opts.species == 'rat'):
        MGImapper = TableEntityMapper(opts.mgi)
        HGNCmapper = TableEntityMapper(opts.hgnc)
        RGDmapper = None
        prefix = 'rgd:'
    else:
        sys.exit("Unrecognized species.  Please choose mouse, rat or human.")


    #read in the query genes and build the request string
    if opts.weighted_query is not None:
        weighted_query=readNodeWeights(opts.weighted_query)
        get_these=weighted_query.keys()
    elif opts.query is not None:
        gene_file = opts.query
        f = open(gene_file, 'r')
        get_these = []
        for line in f:
            get_these.append(line.rstrip())
    else:
        sys.exit("Please provide a query of weighted query file.")

    requestString = " ".join(get_these)

    #establishes connections to ndex.  Search via direct neighbors
    if opts.password is not None and opts.username is not None:
        myNdex = nc.Ndex("http://ndexbio.org", username=opts.username, password=opts.password)
    else:
        myNdex = nc.Ndex("http://ndexbio.org")

    myNet = myNdex.get_neighborhood(opts.network_uuid, requestString, search_depth=1)

    #network=util.ndexPropertyGraphNetworkToNetworkX(myNet)

    #convert Ndex to sif, mapping orthology and collapsing baseterms
    wrapped = NdexToGeneSif(myNet, MGImapper=MGImapper, HGNCmapper=HGNCmapper, RGDmapper=RGDmapper, prefix=prefix)

    wrapped.writeSIF(opts.sif)

    if opts.bel is not None:
        wrapped.writeBELScript(fileName=opts.bel)

    #calculates the heat kernel
    ker = kernel.SciPYKernel(opts.sif, time_T=opts.diffusion_time)

    ker.writeKernel(opts.kernel)

    #establishes and diffuses the query vector
    if opts.query is not None:
        queryVec = queryVector(get_these, ker.labels)
    elif opts.weighted_query is not None:
        queryVec = weightedQueryVector(get_these, ker.labels, weighted_query)

    diffused = ker.diffuse(queryVec)

    sorted_diffused = sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)

    f = open(opts.diffused_query, 'w')
    for key, value in sorted_diffused:
        f.write(str(key)+'\t'+str(value)+'\n')

    f.close()

    #filter the sif, leaving only interactions between the top N genes
    scores = readNodeWeights(opts.diffused_query)
    ints = readSif(opts.sif)
    filtered = filterSif(ints, scores, desired_nodes=opts.entities)

    output = open(opts.filtered_sif, 'w')
    for s in filtered:
        output.write('\t'.join(s) + '\n')


#A=nx.adjacency(network)

#myWrapper = util.NetworkWrapper(myNet)
#myWrapper.writeSummary("RAS-RAF-MEK-MAPK.txt")

#here we are copying util.NetworkWrapper

