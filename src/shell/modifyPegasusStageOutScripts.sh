for dirname in `find ./work/ -type d  -name Downsam\* -print`; do
	echo $dirname;
	for i in `ls $dirname/stage_out*in`;
		do echo $i; foldname=`echo $i|awk -F / '{print $3}'`; echo $foldname; sed 's/\/work\/polyactis\/pegasus\//\/u\/home\/eeskin\/polyacti\/NetworkData\/vervet\/vervetPipeline\/'$foldname'\//'  $i > $i.correctpath ; mv $i.correctpath $i;
	done;
done
