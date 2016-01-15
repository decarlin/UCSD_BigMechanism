__author__ = 'decarlin'

# This script is called from the command line to run the enrichment server with the persisted e_sets
#
# The script reads all of the e_sets and then starts the bottle server
#
# The optional argument 'verbose' specifies verbose logging for testing
#

#
# python run_e_service.py
#
# python run_e_service.py --verbose
#

# body

import uuid
import argparse
from bottle import route, run, template, default_app, request, post, abort
import ndex.client as nc
import json
import kernel_scipy as kernel
import automateHeatKernel as ahk
from multiplyEstablishedKernel import EstablishedKernel
import csv
import operator

parser = argparse.ArgumentParser(description='run the diffusion service')

parser.add_argument('--verbose', dest='verbose', action='store_const',
                    const=True, default=False,
                    help='verbose mode')

arg = parser.parse_args()

if arg.verbose:
    print "Starting diffusion service in verbose mode"
else:
    print "Starting diffusion service"

app = default_app()
app.config['verbose'] = arg.verbose

app.config['ndex'] = nc.Ndex()

current_kernel_id=None
username=None
password=None

@route('/hello')
def index(name='User'):
    verbose_mode = app.config.get("verbose")
    if verbose_mode:
        return template('<b>This is the test method saying Hello verbosely</b>!', name=name)
    else:
        return 'Hello!'

@route('/<network_id>/generate_ndex_heat_kernel', method='GET')
def generate_heat_kernel(network_id):

    if password is not None and ousername is not None:
        myNdex = nc.Ndex("http://ndexbio.org", username=username, password=password)
    else:
        myNdex = nc.Ndex("http://ndexbio.org")
    myNet = myNdex.get_complete_network(network_id)
    wrapped = ahk.NdexToGeneSif(myNet)
    edges=wrapped.edgeList()
    ker = kernel.SciPYKernel(edges)
    kernel_id=uuid.uuid1()
    ker.writeKernel('kernels/{0}'.format(str(kernel_id)))
    
    return json.dumps({'kernel_id':str(kernel_id)})

@route('/generate_network_heat_kernel', metod='POST')
def generate_heat_kernel(network_id):
    dict=json.load(request.body)
    
    ker = kernel.SciPYKernel(dict.get('network'), time_T=opts.diffusion_time)
    

@route('/rank_entities', method='POST')
def diffuse_and_rank():
    dict=json.load(request.body)
    get_these=dict.get('identifier_set')
    kernel_id=dict.get('kernel_id')
    if not get_these:
        abort(401, "requires identifier_set parameter in POST data")
    if not kernel_id:
        abort(401, "requires kernel_id parameter in POST data")

    if kernel_id is not current_kernel_id:
        try:
            ker=EstablishedKernel('kernels/{0}'.format(kernel_id))
        except IOError:
            abort(401, "could not find kernel file matching kernel_id")
    
    queryVec=ahk.queryVector(get_these,ker.labels)
    diffused=ker.diffuse(queryVec)
    sorted_diffused = sorted(diffused.items(), key=operator.itemgetter(1), reverse=True)
    dict_out={"ranked_entities":sorted_diffused}
    return json.dumps(dict_out)

run(app, host='0.0.0.0', port=5602)
