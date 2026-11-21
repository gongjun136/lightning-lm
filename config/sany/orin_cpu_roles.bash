# Source before run_sany_lidar_loc.sh on the inspected 12-core AGX Orin.
# The process starts with the full mask; C++ narrows background threads before
# initialization and explicitly assigns sensor, LIO worker, output, executor and matching.
export LIGHTNING_LM_CPU_AFFINITY=0-11
export LIGHTNING_LM_BACKGROUND_CPUS=0-9
export LIGHTNING_LM_SENSOR_CPUS=10
export LIGHTNING_LM_OUTPUT_CPUS=11
export LIGHTNING_LM_EXECUTOR_CPUS=11
# LIO thread zero stays on CPU10; its second thread shares CPU11 with callbacks/output.
export LIGHTNING_LM_LIO_WORKER_CPUS=${LIGHTNING_LM_LIO_WORKER_CPUS:-11}
export LIGHTNING_LM_LIO_THREADS=${LIGHTNING_LM_LIO_THREADS:-2}
export LIGHTNING_LM_MATCH_CPUS=${LIGHTNING_LM_MATCH_CPUS:-0-9}
export LIGHTNING_LM_NDT_THREADS=${LIGHTNING_LM_NDT_THREADS:-4}
export LIGHTNING_LM_SOLID_CPU_AFFINITY=0-9
# Avoid idle OpenMP teams consuming CPU while the sensor thread waits on them.
export OMP_WAIT_POLICY=PASSIVE
export OMP_PROC_BIND=false
