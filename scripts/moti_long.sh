# delete previous temporary files.
rm /mnt/nvme/indices_upd/bigann/100M_shadow*
rm /mnt/nvme/indices_upd/bigann/100M_merge*
# loop from 1 to 100, merge every 6 iters (18 iters in total), reboot after merge to avoid memory leak.
for i in {0..18}
do
    # Path: build/overall.sh
    # Run the build script
    echo "Running with i = $i"
    OMP_PLACES=cores OMP_PROC_BIND=close build/tests/motivation uint8 /mnt/nvme/data/bigann/bigann.bbin 128 96 1.2 0 0 /mnt/nvme/indices_upd/bigann/100M /mnt/nvme/data/bigann/bigann_query.bbin /mnt/nvme/indices_upd/bigann_gnd/500M_topk 10 4 100 $i 6000000 20 30 40 50 60 80 100
done
