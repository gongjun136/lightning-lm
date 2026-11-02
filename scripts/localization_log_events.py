"""Language-independent event recognition with compatibility for old logs."""
import re

# All strings in each tuple are required for that legacy event. New explicit
# event codes take precedence, so quoted prose cannot change the diagnosis.
LEGACY_EVENTS = {
    'LOC_INITIAL_POSE_UNSET': ('initial pose is not configured',),
    'LOC_MATCH_FAILED': ('localization failed!',),
    'LOC_INIT_MATCH_FAILED': ('初始化匹配失败',),
    'MATCH_CAPTURE_SAVED': ('MATCH_FAILURE_CAPTURE', 'saved='),
    'LOC_OUTPUT_DISABLED': ('localization pose and map outputs disabled',),
    'LOC_OUTPUT_STOPPED': ('localization pose and map outputs stopped',),
    'LOC_OUTPUT_RECOVERED': ('localization pose and map outputs recovered',),
    'LOC_SENSOR_QUEUE_LAG': ('sensor queue lag reached',),
    'LOC_SENSOR_QUEUE_RECOVERED': ('sensor queue lag recovered',),
    'LIO_STALE_CORRECTION_DROPPED': ('drop stale lidar correction before matching:',),
    'LIO_ADAPTIVE_LOAD_CHANGED': ('adaptive lidar load step ',),
    'LOC_RELOCALIZATION_ACCEPTED': ('GLOBAL_RELOCALIZATION[', '] accepted'),
    'SOLID_KISS_REJECT': ('SOLID_KISS_REJECT',),
    'LOC_TRACKING_LOST': ('tracking lost after',),
    'LOC_CAN_SPEED_IGNORED': ('ignore CAN wheel speed:',),
    'IMU_LONG_INTERVAL': ('propagate over long imu interval:',),
    'SOLID_CLOUD_LOAD_FAILED': ('failed to load SOLiD ICP cloud',),
    'SOLID_TARGET_BUILD_FAILED': ('failed to build SOLiD+KISS target',),
    'SOLID_QUERY_PREPARE_FAILED': ('failed to prepare SOLiD+KISS query',),
    'SMOOTHER_MOTION_EXCESSIVE': ('smoother motion is too large:',),
}
EVENT = re.compile(r'(?<!\S)event=([A-Z][A-Z0-9_]*)(?=[\s,]|$)')
CONTEXT_EVENTS = (
    'LOC_INITIAL_POSE_UNSET', 'LOC_MATCH_FAILED', 'LOC_INIT_MATCH_FAILED', 'MATCH_CAPTURE_SAVED', 'LOC_OUTPUT_DISABLED',
    'LOC_OUTPUT_STOPPED', 'LOC_OUTPUT_RECOVERED', 'LOC_SENSOR_QUEUE_LAG',
    'LOC_SENSOR_QUEUE_RECOVERED', 'LIO_STALE_CORRECTION_DROPPED',
    'LOC_CAN_SPEED_IGNORED', 'IMU_LONG_INTERVAL', 'SOLID_CLOUD_LOAD_FAILED',
    'SOLID_TARGET_BUILD_FAILED', 'SOLID_QUERY_PREPARE_FAILED',
)
BOUNDARY_EVENTS = (
    'LOC_MATCH_FAILED', 'LOC_INIT_MATCH_FAILED', 'MATCH_CAPTURE_SAVED', 'LOC_OUTPUT_DISABLED', 'LOC_OUTPUT_RECOVERED',
    'LOC_INITIAL_POSE_UNSET', 'LOC_SENSOR_QUEUE_LAG', 'LOC_RELOCALIZATION_ACCEPTED',
)


def has_event(text, code):
    explicit = EVENT.search(text)
    if explicit:
        return explicit[1] == code
    return all(marker in text for marker in LEGACY_EVENTS[code])


def any_event(text, codes):
    explicit = EVENT.search(text)
    if explicit:
        return explicit[1] in codes
    return any(all(marker in text for marker in LEGACY_EVENTS[code]) for code in codes)
