#!/usr/bin/python3

import shlex
from subprocess import Popen, PIPE,call,check_output
import argparse,datetime,subprocess


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description="Build Docker file for conditional analysis")

    parser.add_argument("--image", type= str,
                        help="name of image",default = 'conditional_analysis')
    parser.add_argument("--version", type= str,
                        help="version value, e.g.0.001",required = True)
    parser.add_argument("--push",action = 'store_true')
    parser.add_argument("--args",type = str,default = '')
    args = parser.parse_args()

    
    basic_cmd = 'docker build -t eu.gcr.io/finngen-refinery-dev/' + args.image +':' +args.version
    cmd = basic_cmd + ' -f Dockerfile ..' + ' ' + args.args
    print(cmd)
    call(shlex.split(cmd))

    if args.push:
        cmd = ' docker -- push eu.gcr.io/finngen-refinery-dev/' + args.image +':' + args.version
        print(cmd)
        call(shlex.split(cmd))
