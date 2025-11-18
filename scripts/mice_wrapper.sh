k=31
filelist=TODO
workdir=test
ggcat build -k $k -s 1 -l $filelist -o $workdir/graph.gfa --gfa-v1
gfa2gff $workdir/graph.gfa $(cat $filelist)  > $workdir/graph.gff
#rm $workdir/graph.gfa
mice $workdir/graph.gff -o $workdir/mice_out -r 2
#rm $workdir/graph.gff
