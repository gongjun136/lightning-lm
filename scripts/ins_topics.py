#!/usr/bin/env python3
"""Print CGI input topics, or mode-specific outputs for recording."""
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

def outputs(filename):
    c=yaml.safe_load(pathlib.Path(filename).read_text())
    mode=c['ins_only'].get('output_reference','device_output')
    if mode=='device_output':
        return ['/localization/ins_reference/'+x for x in
                ['pose','status','fault_status','path','diagnostics']]
    if mode!='rear_axle':
        raise ValueError('未知 INS 输出参考点，请声明 device_output 或 rear_axle')
    return ['/slamPoseRaw_topic','/localization/pose_vel','/localization/loc_status',
            '/localization/fault_status','/localization/ins_diagnostics',
            '/diagnostics/heartbeat/lightning_slam']

if __name__=='__main__':
    print('\n'.join(outputs(sys.argv[2]) if sys.argv[1]=='--outputs' else topics(sys.argv[1])))
