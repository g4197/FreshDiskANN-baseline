# bs = 30M
for i in {0..7}
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PLACES=cores OMP_PROC_BIND=close build/tests/motivation_stress uint8 /mnt/nvme/data/bigann/bigann.bbin 128 128 1.2 0 0 /mnt/nvme3/indices/bigann/800M /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd_insert/1B_topk 10 4 200 $i 30000000 25
    # build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann_2M.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/1m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/2M_topk 10 4 100 $i 20000 20 30 40
done
