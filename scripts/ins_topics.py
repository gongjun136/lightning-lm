#!/usr/bin/env python3
"""Print the exact input topic allowlist for CGI-only recording/replay."""
import pathlib
import sys
import yaml

def topics(filename):
    filename=pathlib.Path(filename)
    c=yaml.safe_load(filename.read_text())
    prefix=c['ins_only'].get('topic_prefix','/cgi430')
    result=[prefix+'/'+x for x in ['time','position/latitude','position/longitude',
        'position/altitude','position/sigma','attitude','attitude/sigma',
        'velocity','velocity/sigma','ins/status','driver_status']]
    cloud=c.get('cloud',{})
    if cloud.get('enabled',False):
        base=c
        if 'lidar_config' in cloud:
            base=yaml.safe_load((filename.parent/cloud['lidar_config']).read_text())
        result += [v for k,v in base['multi_lidar']['topics'].items() if k.startswith('lidar_')]
    return result

if __name__=='__main__':
    print('\n'.join(topics(sys.argv[1])))
