# loop from 1 to 100
for i in {0..18}
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PROC_BIND=close build/tests/overall_performance uint8 /mnt/nvme/data/bigann/bigann.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/100m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/500M_topk 10 4 100 $i 6000000 20 30 32 35 37 40 45
done
