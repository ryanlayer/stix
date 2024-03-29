import subprocess
import json
import sys
from flask import Flask,request,render_template
import toolshed as ts
from operator import itemgetter
import collections
import argparse

app = Flask(__name__)
stix_path = None
index_path = None
db_path = None

@app.route('/')
def hello():
    raw = request.args.get('raw')
    left_bp = request.args.get('left')
    right_bp = request.args.get('right')
    region =  request.args.get('region')
    sv_type = request.args.get('sv_type')


    if sv_type != None:
        if region != None:
            try:
                chrm = region.split(':')[0]
                left = int(region.split(':')[1].split('-')[0])
                right = int(region.split(':')[1].split('-')[1])
            except:
                return '{"error":"unsupported variant format"}'

            left_bp = chrm + ':' + str(left) + '-' + str(left + 1)
            right_bp = chrm + ':' + str(right) + '-' + str(right + 1)
        elif left_bp == None or right_bp == None:
            return '{"error":"unsupported variant format"}'

        print(left_bp, right_bp)

    if raw != None:
        print(stix_path, sv_type, index_path, db_path)
        cmd = stix_path + \
            ' -t ' + sv_type + \
            ' -i ' + index_path + \
            ' -d ' + db_path + \
            ' -s 500 -j '

        print cmd

        if ((left_bp != None) and (right_bp != None)):
            cmd += ' -l ' + left_bp
            cmd += ' -r ' + right_bp
            print request.remote_addr, cmd
            proc = subprocess.Popen(cmd,
                                    shell=True,
                                    stderr=sys.stderr,
                                    stdout=subprocess.PIPE)
            output = proc.stdout.read()
            proc.wait()

            return output
        else:
            return "{}"
    else:
        sv = { 'sv_type': sv_type,
               'left' : left_bp,
               'right' : right_bp }
        return render_template('sv_view.html', data=sv)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='STIX web server.')

    parser.add_argument('--stix',
                        dest='stix_path',
                        help='Path to stix')

    parser.add_argument('--port',
                        dest='port',
                        help='Port to run on')

    parser.add_argument('--db',
                        dest='db',
                        help='STIX ped db')

    parser.add_argument('--index',
                        dest='index',
                        help='STIX index')


    args = parser.parse_args()

    stix_path = args.stix_path
    port = args.port
    index_path = args.index
    db_path = args.db
    print(db_path)
    app.run(host='0.0.0.0', port=port)
