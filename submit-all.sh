# this is a test:
for START_INDEX in 0
do
    echo $START_INDEX
    export START_INDEX
    qsub -V -qshortq run.qsub
done

# for START_INDEX in {0..99}
# do
#     echo $START_INDEX
#     export START_INDEX
#     qsub -V -qshortq run.qsub
# done
