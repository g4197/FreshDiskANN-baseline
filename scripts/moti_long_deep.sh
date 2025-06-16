# loop from 1 to 100
for i in {0..18}
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PLACES=cores OMP_PROC_BIND=close build/tests/motivation float /mnt/nvme/data/deep/200M.fbin 128 64 1.2 0 0 /mnt/nvme/indices/deep/100M /mnt/nvme/data/deep/queries.fbin /mnt/nvme/data/deep/gnd_insert/200M_topk 10 4 100 $i 6000000 25 30 33 37 40 45
    # build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann_2M.bbin 128 64 1.2 0 0 /mnt/nvme/indices/bigann/1m /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/data/bigann/gnd/2M_topk 10 4 100 $i 20000 20 30 40
done
