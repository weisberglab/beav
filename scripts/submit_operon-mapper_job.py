#!/usr/bin/env python
import requests
import re
import getopt, sys

#email="weisbeal@oregonstate.edu"
email=""
strain=""
#strain="AS1D4"


try:
    opts, args = getopt.getopt(sys.argv[1:], "es", ["email=", "strain="])
except getopt.GetoptError as err:
    # print help information and exit:
    print(err)  # will print something like "option -a not recognized"
    sys.exit(2)
for o, a in opts:
    if o in ("-e", "--email"):
        email = a
    elif o in ("-s", "--strain"):
        strain = a
    else:
        assert False, "unhandled option"

#in case strain name is a number
strain=str(strain)

userdata = {"email1": email, "descri": strain + " operons", "genepairs": "si", "operons": "si", "cogs": "si","orfsdescri": "si"}
files = {'fastafile': open(strain + '.fna', 'rb'), 'gfffile': open(strain + '.gbff.gff', 'rb')}

try:
	#need form action url for post, not original webpage
    resp = requests.post('https://biocomputo.ibt.unam.mx/operon_mapper/capta_forma_01.pl', data=userdata, files=files, headers={'referer': 'https://biocomputo.ibt.unam.mx/', 'origin': 'https://biocomputo.ibt.unam.mx'})
    resp.raise_for_status()
except requests.exceptions.HTTPError as errh:
    print ("Http Error:",errh)
    sys.exit(0)
except requests.exceptions.ConnectionError as errc:
    print ("Error Connecting:",errc)
    sys.exit(0)
except requests.exceptions.Timeout as errt:
    print ("Timeout Error:",errt)
    sys.exit(0)
except requests.exceptions.RequestException as err:
    print ("Error:",err)
    sys.exit(0)

with open("testresponse.html", "w") as outfile:
    outfile.write(resp.text)
outfile.close()

with open("operon-mapper_results_url", "w") as resfile:
    resfile.write(str(resp.status_code) + '\n')

    for line in resp.text.splitlines():
        if '<meta http-equiv="refresh"' in line: 
            new_urltemp=re.sub(r'^.*url=', '', line)
            new_url=re.sub(r'".*', '', new_urltemp)
            resfile.write(str(new_url) + '\n')
resfile.close()

#<input type="checkbox" name="genepairs" value="si" >
#Predicted operonic gene pairs

#<input type="checkbox" name="operons" value="si" >
#Predicted operons

#<input type="checkbox" name="orfcoordinates" value="si" >
#Predicted ORFs coordinates

#<input type="checkbox" name="dna" value="si" >
#DNA sequences of the predicted ORFs

#<input type="checkbox" name="protein" value="si" >
#Protein sequences of the translated predicted ORFs

#<input type="checkbox" name="cogs" value="si" >
#COGs assignations

#<input type="checkbox" name="orfsdescri" value="si" >

#<input type="checkbox" name="todos" value="si" >
#All possible outfiles

#<input type="checkbox" name="todoscomprimido" value="si"  checked>
#All possible outfiles and a compressed file with all of them
