#!/usr/bin/env python

__author__ = 'danielcarlin'


import automateHeatKernel as ahk
import operator
from optparse import OptionParser
import sys
import ndex.client as nc
import kernel_scipy as kernel
import csv
import weighted_kernel as wk


if __name__ == "__main__":
    parser = OptionParser()

    parser.add_option("-k", "--kernel", dest="kernel", action="store", type="string", default="kernel.txt",
                      help="Output file for kernel")
    parser.add_option("-t", "--diffusion-time", dest="diffusion_time", type="float", default=0.1,
                      help="Time parameter (scale) for diffusion.  Default is 0.1")
    parser.add_option("-d", "--diffused-query", dest="diffused_query", type="string", default="diffused.txt",
                      help="Output file for diffused query.")
    parser.add_option("-u", "--username", dest="username", type="string", default=None, help="NDEx username")
    parser.add_option("-p", "--password", dest="password", type="string", default=None, help="NDEx password")
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
    parser.add_option("-n", "--network-uuid", dest="network_uuid", type="string",
                      default="9ea3c170-01ad-11e5-ac0f-000c29cb28fb",
                      help="uuid of query, default is for the BEL large corpus")
    parser.add_option("-s", "--sif", dest="sif", action="store", type="string", default="out.sif",
                      help="Output file for intermediate sif, translated to HGNC")


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


    #establishes connections to ndex.  Search via direct neighbors
    if opts.password is not None and opts.username is not None:
        myNdex = nc.Ndex("http://ndexbio.org", username=opts.username, password=opts.password)
    else:
        myNdex = nc.Ndex("http://ndexbio.org")

    myNet = myNdex.get_complete_network(opts.network_uuid)

    wrapped = ahk.NdexToGeneSif(myNet, MGImapper=MGImapper, HGNCmapper=HGNCmapper, RGDmapper=RGDmapper, prefix=prefix)

    wrapped.writeSIF(opts.sif)

    ker = kernel.SciPYKernel(opts.sif, time_T=opts.diffusion_time)

    ker.writeKernel(opts.kernel)

