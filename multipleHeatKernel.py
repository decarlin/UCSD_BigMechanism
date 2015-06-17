#!/usr/bin/env python

import automateHeatKernel as ahk
import operator
from optparse import OptionParser
import sys
import ndexClient as nc
import kernel_scipy as kernel
import csv
import weighted_kernel as wk


if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option("-q", "--query", dest="query", action="store", type="string",
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
    parser.add_option("-n", "--networks-file", dest="network_uuids_file", type="string",
                      default="network_uuids.txt",
                      help="uuids of query networks")
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
    parser.add_option("-c", "--combination-method", dest="comb", type="string", default="sif", help="Method of combining different networks.  Default is 'sif', which creates a single .sif of all networks in the target namespace. Another option is 'kernel', which creates a heat kernel for each network, then multiplies by each in turn, which has the effect of double counting shared edges."  )
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
    if (opts.species == 'human'):
        MGImapper = ahk.TableEntityMapper(opts.mgi)
        RGDmapper = ahk.TableEntityMapper(opts.rgd)
        HGNCmapper = None
        prefix = 'hgnc:'
    elif (opts.species == 'mouse'):
        HGNCmapper = ahk.TableEntityMapper(opts.hgnc)
        RGDmapper = ahk.TableEntityMapper(opts.rgd)
        MGImapper = None
        prefix = 'mgi:'
    elif (opts.species == 'rat'):
        MGImapper = ahk.TableEntityMapper(opts.mgi)
        HGNCmapper = ahk.TableEntityMapper(opts.hgnc)
        RGDmapper = None
        prefix = 'rgd:'
    else:
        sys.exit("Unrecognized species.  Please choose mouse, rat or human.")

    gene_file = opts.query

    f = open(gene_file, 'r')

    get_these = []

    for line in f:
        get_these.append(line.rstrip())

    requestString = " ".join(get_these)

    #establishes connections to ndex.  Search via direct neighbors
    if opts.password is not None and opts.username is not None:
        myNdex = nc.Ndex("http://ndexbio.org", username=opts.password, password=opts.username)
    else:
        myNdex = nc.Ndex("http://ndexbio.org")

    fn = open(opts.network_uuids_file)

    if opts.comb == 'sif':
        num_networks=0.0;
        for net_line in fn:
            myNet = myNdex.getNeighborhood(net_line.rstrip(), requestString, searchDepth=1)

            wrapped = ahk.NdexToGeneSif(myNet, MGImapper=MGImapper, HGNCmapper=HGNCmapper, RGDmapper=RGDmapper, prefix=prefix)

            wrapped.writeSIF(opts.sif, append=True)
            num_networks=num_networks+1

        ker = kernel.SciPYKernel(opts.sif, time_T=opts.diffusion_time)

        ker.writeKernel(opts.kernel)

        #establishes and diffuses the query vector
        queryVec = ahk.queryVector(get_these, ker.labels)

        diffused = ker.diffuse(queryVec)
        sorted_diffused = sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)
    elif opts.comb =='kernel':
        num_networks=0.0;
        for net_line in fn:
            myNet = myNdex.getNeighborhood(net_line.rstrip(), requestString, searchDepth=1)

            wrapped = ahk.NdexToGeneSif(myNet, MGImapper=MGImapper, HGNCmapper=HGNCmapper, RGDmapper=RGDmapper, prefix=prefix)

            wrapped.writeSIF('multi.sif', append=True)

            num_networks=num_networks+1

        wk.MultiSifToWeightedSif('multi.sif',transform=wk.log2plus1, outfile=opts.sif)

        ker = wk.WeightedSciPYKernel(opts.sif, time_T=(opts.diffusion_time/num_networks))

        ker.writeKernel(opts.kernel)
        #establishes and diffuses the query vector
        queryVec = ahk.queryVector(get_these, ker.labels)

        diffused = ker.diffuse(queryVec)
        sorted_diffused = sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)

    writer = csv.writer(open(opts.diffused_query, 'wb'), delimiter='\t')
    for key, value in sorted_diffused:
        writer.writerow([key, value])

        #filter the sif, leaving only interactions between the top N genes
    scores = ahk.readNodeWeights(opts.diffused_query)
    ints = ahk.readSif(opts.sif)
    filtered = ahk.filterSif(ints, scores, desired_nodes=opts.entities)

    output = open(opts.filtered_sif, 'wb')
    for s in filtered:
        output.write('\t'.join(s) + '\n')
