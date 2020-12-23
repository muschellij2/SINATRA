SHAPES_PER_CLASS=(25 50 75)
NUM_CONES=(25 50 75)
DIR_PER_CONE=(5 10)
CURVE_LENGTH=(25 50)

for i in {0..2}
 do
 	for j in {0..2}
 	do
 		for k in {0,1}
 		do
 			for l in {0,1}
 			do
				task_id=$(( $(( 2 * $(( $(( 2 * $(( $(( 3 * $(( $i )) )) + $j )) )) + $k )) )) + $l ))
				shapes=${SHAPES_PER_CLASS[$i]}
				numcones=${SHARED[$j]}
				dirpercone=${DIR_PER_CONE[$k]}
				curvelength=${CURVE_LENGTH[$l]}


				echo $task_id
			done
		done
	done
done