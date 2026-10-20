# Source before run_sany_lidar_loc.sh on the inspected 12-core AGX Orin.
# The process starts with the full mask; C++ narrows background threads before
# initialization and explicitly assigns only the three critical role threads.
export LIGHTNING_LM_CPU_AFFINITY=0-11
export LIGHTNING_LM_BACKGROUND_CPUS=0-9
export LIGHTNING_LM_SENSOR_CPUS=10
export LIGHTNING_LM_OUTPUT_CPUS=11
export LIGHTNING_LM_EXECUTOR_CPUS=11
export LIGHTNING_LM_SOLID_CPU_AFFINITY=0-9
# Avoid idle OpenMP teams consuming CPU while the sensor thread waits on them.
export OMP_WAIT_POLICY=PASSIVE
export OMP_PROC_BIND=false
